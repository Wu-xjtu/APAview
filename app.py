"""APAview Flask app (server).
Usage:
	# Initialize Flask Python environment yourself here
	cd /path/of/this/file
	# For Windows Powershell
	$env:FLASK_ENV="development"
	# Or for Linux bash
	# export FLASK_ENV=development
	flask run
"""

### Usage agreement
# If you are interested in obtaining the software for commercial use, please contact Xiaoming Wu (wxm@mail.xjtu.edu.cn). For academic use, downloading or using the software means you agree with the MIT License.
### Citation
# Xi Hu, Jialin Song, Jinping Wan, Jianqiang Du, Junbo Duan, Huqin Zhang, Xiaoming Wu, APAview: a web-based platform for alternative polyadenylation analyses in hematological tumors


import os
import math
import operator
from collections import namedtuple
import re
import logging

# Third-party imports
from flask import Flask, render_template, request, jsonify, make_response, abort, redirect
import flask.logging
import sqlite3
from werkzeug.exceptions import BadRequest

import numpy as np

from pandas import Series, DataFrame
import pandas as pd

import xarray as xr

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
plotly_colors = px.colors.qualitative.Plotly

import statsmodels.api as sm
from statsmodels.stats.weightstats import ttest_ind
from statsmodels.stats.oneway import anova_oneway

# Debug imports
# import pickle
# import pdb

# Global vars (will be loaded in app.before_first_request)
db = None
APA_datasets = {}	# Containing all datasets

# Global constants (database column names and others)
config = {
	"sql_file":"APA.sqlite",
	"gene_info_sql_file":"gene_info.sqlite",
	"data_file":"data.nc",
	"max_items":499,
	"ref_genome":{"QAPA":"hg38","DaPars":"hg19"},
	# "valid_case_min_ratio":0.5,
	"valid_case_min":3,
	# "case_show_columns":["rowid","sample_id","gender","age","diagnosis","is_recurrence","duplication_rate","cancer_cell_rate","fusion_genes","sample_source"],
	"diseases":{"AA","AL","ALL","AML","APL","B-ALL","MDS","MPAL","normal","T-ALL","T-LBL"},
	"default_disease":"AML",
	"QAPA_show_columns":["id","gene_id","gene_symbol","strand","chrom","expressed_case_count",
		("{}||' - '||{} AS last_exon_region", ("last_exon_start", "last_exon_end")),
		("{}||' - '||{} AS UTR3_region", ("UTR3_start", "UTR3_end")),
		("{}-{} AS length", ("UTR3_end", "UTR3_start")),
		("CASE WHEN {}='+' THEN {} ELSE {} END AS UTR3_terminal", ("strand", "UTR3_end", "UTR3_start"))
	],
	"DaPars_show_columns":["id","gene_id","gene_symbol","strand","chrom","expressed_case_count","{0}_truncation_site AS truncation_site",
		("{}||' - '||{} AS UTR3_region", ("UTR3_start", "UTR3_end")),
		("{}-{} AS length", ("UTR3_end", "UTR3_start")),
		("CASE WHEN {}='+' THEN {} ELSE {} END AS UTR3_terminal", ("strand", "UTR3_end", "UTR3_start")),
		("CASE WHEN {}='+' THEN {}-{} ELSE {}-{} END AS truncation_length", ("strand", "{0}_truncation_site", "UTR3_start", "UTR3_end", "{0}_truncation_site"))
	]
}


# Create global Flask app (see Flask documentation)
flask.logging.default_handler.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)s in %(funcName)s(): %(message)s"))
app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 4 * 2**20


@app.before_first_request
def initialize():
	"""Connect to database and datasets (initialize global variables)."""
	global db

	app.logger.info("Connecting to database")
	db = sqlite3.connect(config["sql_file"], check_same_thread=False)
	db.set_trace_callback(app.logger.debug)
	db.row_factory = sqlite3.Row	# Return list of dict, not list of tuple
	db.execute("ATTACH DATABASE ? AS gene_info_db", (config["gene_info_sql_file"],))
	db.execute("PRAGMA case_sensitive_like=true")
	db.execute("PRAGMA query_only=true")
	config["n_cases"] = {}
	config["n_cases"]["QAPA"] = db.execute("SELECT max(rowid) FROM QAPA_case_info").fetchall()[0][0]
	config["n_cases"]["DaPars"] = db.execute("SELECT max(rowid) FROM DaPars_case_info").fetchall()[0][0]
	app.logger.info("Successfully connecting to database.")
	app.logger.info("Start initializing HDF5 data")
	APADataset = namedtuple("APADataset", ["PDUI", "TPM", "site_TPM"])	# Pack PDUI and TPM dataset from the same source into this
	data_file_name = config["data_file"]
	# TODO: Is it necessary to open the data file multiple times separately?
	APA_datasets["case_type"] = xr.open_dataarray(data_file_name,group="/case_type",engine="h5netcdf")
	# xarray DataArrays are lazily loaded from disk when indexed. This is what we need.
	PDUI_dataset = xr.open_dataarray(data_file_name,group="/QAPA/PDUI",engine="h5netcdf")
	TPM_dataset = xr.open_dataarray(data_file_name,group="/QAPA/logTPM_gene",engine="h5netcdf")
	site_TPM_dataset = xr.open_dataarray(data_file_name,group="/QAPA/logTPM_site",engine="h5netcdf")
	APA_datasets["QAPA"] = APADataset(PDUI=PDUI_dataset, TPM=TPM_dataset, site_TPM=site_TPM_dataset)
	# For DaPars data, there's only one APA site on each gene.
	PDUI_dataset = xr.open_dataarray(data_file_name,group="/DaPars/PDUI",engine="h5netcdf")
	TPM_dataset = xr.open_dataarray(data_file_name,group="/DaPars/logTPM_gene",engine="h5netcdf")
	APA_datasets["DaPars"] = APADataset(PDUI=PDUI_dataset, TPM=TPM_dataset, site_TPM=TPM_dataset)
	app.logger.info("HDF5 data initializing succeed. Initialization finished.")


# Routing of web pages
@app.route("/", methods=["GET"])
def root():
	return redirect("/index")

@app.route("/index", methods=["GET"])
def index():
	return render_template("index.html.j2")

@app.route("/help", methods=["GET"])
def help():
	return render_template("help.html.j2")

@app.route("/contact", methods=["GET"])
def contact():
	return render_template("contact.html.j2")

@app.route("/<any(QAPA,DaPars):dataset>/query", methods=["GET"])
def query(dataset):
	# FIXME: Is it necessary to treat URL as case-insensitive?
	return render_template("query.html.j2", dataset=dataset, config=config)

@app.route("/de_analysis", methods=["GET"])
def DE_analysis():
	return render_template("de_analysis.html.j2")


@app.after_request
def postprocess_bgzip_range(response):
	"""For a patch when IGV.js requesting .bed.gz data segments.
		"Content-Encoding: gzip" will let the browser automatically decompress the segment and lead to an error.
	"""
	if response.content_range and response.content_encoding == "gzip":
		del response.content_encoding
	return response


def make_json_response(o):
	if not isinstance(o, str):
		o = jsonify(o)
	response = make_response(o)
	response.headers["Content-Type"] = "application/json"
	return response

@app.route("/<any(QAPA,DaPars):dataset>/de_analysis_table", methods=["GET"])
def DE_analysis_table(dataset):
	"""(This URL is invisible to user)
		Return case info table as JSON to bootstrap-table.
	"""
	# data = db.execute("SELECT {} FROM case_info".format(",".join(config["case_show_columns"]))).fetchall()
	data = db.execute("SELECT rowid,* FROM {}_case_info ORDER BY rowid".format(dataset)).fetchall()
	return make_json_response([dict(row) for row in data])


@app.route("/<any(QAPA,DaPars):dataset>/de_analysis_table_header", methods=["GET"])
def DE_analysis_table_header(dataset):
	"""(This URL is invisible to user)
		Return case info table header (columns) as JSON to bootstrap-table.
	"""
	info = db.execute("PRAGMA table_info({}_case_info)".format(dataset)).fetchall()
	column_info = [{"checkbox":True}, {"field":"group", "title":"Group", "visible":True, "formatter":"groupNameFormatter"}, {"field":"rowid", "title":"No.", "visible":True}]
	column_info.extend({"field":l[1], "title":l[1], "visible":l[1] in ("rowid", "sample_id", "gender", "age", "diagnosis")} for l in info)
	return make_json_response(column_info)


def identifier_type(s, raise_error=True):
	"""Determine the type of a gene symbol (identifier)."""
	if re.fullmatch(r"ENS[A-Z]\d{11}(?:\.\d+)?", s):
		return "ensembl"
	elif re.fullmatch(r"ENS[A-Z]\d{11}_\d+_[SPD]", s):
		return "ensembl_apa"
	elif re.fullmatch(r"[NX][A-Z]_\d+(?:\.\d+)?", s):
		return "refseq"
	elif s.isdigit():
		return "entrez"
	elif re.fullmatch(r"[A-Z0-9\-]+", s):
		return "symbol"
	elif raise_error:
		abort(400, description="Can't infer the ID type: \"{}\"".format(s))
	else:
		return "unknown"

def gene_set_identifier_type(s, raise_error=True):
	if re.fullmatch(r"GO:\d{7}", s):
		return "go"
	elif re.fullmatch(r"(?:map|hsa)\d{5}", s):
		return "kegg"
	elif raise_error:
		abort(400, description="Can't infer the gene set ID type: \"{}\"".format(s))
	else:
		return "unknown"

@app.route("/<any(QAPA,DaPars):dataset>/query_table", methods=["GET", "POST"])
def query_table(dataset):
	"""(This URL is invisible to user)
		Query and return APA site table as JSON to bootstrap-table.
		Query parameters are encoded in URL as key-value pairs.
	"""
	# Query parameters (GET) example: 
	# ?IDType=GENE&ID=ENSG00000000001&siteType=SINGLE&chromosome=chr3&regionType=UTR&regionStart=12&regionStop=1200&orderBy=DEFAULT&descendOrder=true
	# Or a POST request containing a ID list file
	return query_table_impl(dataset)


@app.route("/<any(QAPA,DaPars):dataset>/query_table_simple", methods=["GET"])
def query_table_simple(dataset):
	"""(This URL is invisible to user)
		Query and return APA site table as JSON to bootstrap-table.
		Query parameters are encoded in URL as key-value pairs.
	"""
	# Query parameters (GET) example: 
	# ?searchString=ENSG00000000001
	# Detect meaning of query string
	query_str = request.args.get("searchString", "")
	if not query_str:
		return query_table_impl(dataset, {})
	if gene_set_identifier_type(query_str, raise_error=False) != "unknown":
		# Gene set
		return query_table_impl(dataset, {"geneSetID":query_str})
	is_location = re.fullmatch(r"([A-Za-z0-9.]+):(\d*)-(\d*)", query_str)
	if is_location:
		return query_table_impl(dataset, {"chromosome":is_location.group(1), "regionStart":is_location.group(2), "regionStop":is_location.group(3)})
	id_type = identifier_type(query_str.split(";",1)[0], raise_error=False)
	if id_type != "unknown":
		# Gene or APA ID
		is_APA = (dataset == "QAPA" and id_type == "ensembl_apa"
			or dataset == "DaPars" and id_type == "refseq")
		return query_table_impl(dataset, {"ID":query_str, "IDType":"APA" if is_APA else "GENE"})
	else:
		abort(400, description="Can't get the meaning of query string: \"{}\"".format(query_str))


def remove_version(ID):
	return ID.split(".", 1)[0]

def request_assert(condition, description=None):
	if not condition:
		abort(400, description=description)

def query_table_impl(dataset, args=None):
	# "request" is also passed in
	if request.method == "POST":
		file = request.files["IDListFile"]
		ID = [line.strip() for line in file.read().decode("ascii").split("\n")]
		ID = [line for line in ID if line and not line.startswith("#")]
		args = request.form
	else:
		if args is None:
			args = request.args
		ID = args.get("ID", "")
		if ";" in ID:
			ID = ID.split(";")
			ID = [line for line in ID if line]
	if isinstance(ID, str):
		ID = remove_version(ID)
	else:
		ID = [remove_version(s) for s in ID]
		if len(ID) == 1:
			ID = ID[0]

	select_statement = "SELECT {} FROM {} AS d".format(
		",".join("d."+col if isinstance(col, str) else col[0].format(*("d."+s for s in col[1])) for col in config[dataset+"_show_columns"]), 
		dataset+"_event")

	# Choose truncation site type (currently it's a hack, see DaPars_show_columns)
	if dataset == "DaPars":
		disease = args.get("disease", "")
		if not disease:
			disease = config["default_disease"]
		request_assert(disease in config["diseases"])
		select_statement = select_statement.format(disease.replace("-","_"))

	where_statements = []
	params = []

	ID_type = args.get("IDType", "")
	ID_type_enum = {
		"":"id",
		"APA":"id",
		"GENE":"gene_id"
	}

	join_parameters = set()
	join_statement_template = "INNER JOIN {} ON {}={}"
	need_ID_type = "ensembl"
	if ID:	# Not empty string
		request_assert(ID_type in ID_type_enum)
		ID_type = "d."+ID_type_enum[ID_type]

		if ID_type == "d.gene_id":
			if isinstance(ID, str):
				gene_ID_type = identifier_type(ID)
			else:
				gene_ID_type = identifier_type(ID[0])
				request_assert(all(identifier_type(s) == gene_ID_type for s in ID[1:]), "Type of all gene IDs should be the same.")
			app.logger.debug("Querying by ID: inferred identifier type: {}".format(gene_ID_type))
			if not gene_ID_type == need_ID_type:
				app.logger.debug("Identifier convertion is used")
				join_parameters.add(("id2{}".format(need_ID_type), "d.gene_id", "id2{}.accession".format(need_ID_type)))
				join_parameters.add(("id2{}".format(gene_ID_type), "id2{}._id".format(need_ID_type), "id2{}._id".format(gene_ID_type)))
				ID_type = "id2{}.accession".format(gene_ID_type)
		if isinstance(ID, str):
			where_statements.append("{}=?".format(ID_type))
			params.append(ID)
		else:
			where_statements.append("{} IN ({})".format(ID_type, ",".join(["?"]*len(ID))))
			params.extend(ID)

	gene_set = args.get("geneSetID", "")
	if gene_set:
		gene_set_type = gene_set_identifier_type(gene_set)
		join_parameters.add(("id2{}".format(need_ID_type), "d.gene_id", "id2{}.accession".format(need_ID_type)))
		join_parameters.add(("id2{}".format(gene_set_type), "id2{}._id".format(need_ID_type), "id2{}._id".format(gene_set_type)))
		if gene_set_type == "kegg":
			# Get the digital part only
			gene_set = gene_set[-5:]
			where_statements.append("id2{}.gene_set=?".format(gene_set_type))
			params.append(gene_set)
		else:
			assert gene_set_type == "go"
			params.append(gene_set)
			join_parameters.add(("go_term AS go_term1", "id2go.gene_set", "go_term1.go_accession"))
			join_parameters.add(("go_offspring", "go_offspring._offspring_go_id", "go_term1._go_id"))
			join_parameters.add(("go_term AS go_term2", "go_offspring._go_id", "go_term2._go_id"))
			where_statements.append("go_term2.go_accession=?")

	join_statement = " ".join([join_statement_template.format(*p) for p in join_parameters])

	site_type = args.get("siteType", "")
	site_type_enum = {
		"":"",
		"PROXIMAL":"P",
		"DISTAL":"D",
		"SINGLE":"S"
	}
	request_assert(site_type in site_type_enum)
	site_type = site_type_enum[site_type]
	if site_type:
		where_statements.append("d.event_type=?")
		params.append(site_type)

	region_type = args.get("regionType", "")
	region = (args.get("regionStart", ""), args.get("regionStop", ""))
	if region_type == "UTR_TERMINAL" or not region_type:
		start_col_name = end_col_name = "CASE WHEN d.strand='+' THEN d.UTR3_end ELSE d.UTR3_start END"
	else:
		region_type_enum = {
			"UTR":"UTR3",
			"LAST_EXON":"last_exon"
		}
		request_assert(region_type in region_type_enum)
		region_type = region_type_enum[region_type]
		start_col_name = "d."+region_type+"_start"
		end_col_name = "d."+region_type+"_end"
	chromosome = args.get("chromosome", "")
	if chromosome and not chromosome.startswith("chr"):
		chromosome = "chr"+chromosome

	if chromosome or region[0] or region[1]:	# Chromosome must be non-empty when regions are non-empty
		where_statements.append("d.chrom=?")
		params.append(chromosome)
	if region[0] and region[1]:	# non-empty
		where_statements.append("{} BETWEEN ? AND ?".format(start_col_name))
		params.extend((int(region[0]), int(region[1])))
		if end_col_name != start_col_name:
			last_i = len(params)
			where_statements.append("{} BETWEEN ?{} AND ?{}".format(end_col_name, last_i-1, last_i))
	elif region[0]:
		where_statements.append("{}>=?".format(start_col_name))
		params.append(int(region[0]))
	elif region[1]:
		where_statements.append("{}<=?".format(end_col_name))
		params.append(int(region[1]))

	if where_statements:	# Not empty
		where_statement = "WHERE " + " AND ".join(where_statements)
	else:
		where_statement = ""

	orderby = args.get("orderBy", "")
	desc_order = args.get("descendOrder", "")
	request_assert(desc_order in ("", "true", "false"))
	desc_order = desc_order == "true"
	order_statement = "DESC" if desc_order else "ASC"
	orderby_enum = {
		"":"",
		"DEFAULT":"",
		"GENE_AND_NUMBER":["d.gene_id", "d.number"],
		"LOCATION":["d.chrom", end_col_name if desc_order else start_col_name]
	}
	request_assert(orderby in orderby_enum)
	orderby = orderby_enum[orderby]
	if orderby:
		orderby_statement = "ORDER BY "+", ".join("{} {}".format(s, order_statement) for s in orderby)
	else:
		orderby_statement = ""

	limit_statement = "LIMIT ?"
	params.append(config["max_items"]+1)
	sql_statement = " ".join([select_statement, join_statement, where_statement, orderby_statement, limit_statement])

	data = db.execute(sql_statement, params).fetchall()
	return make_json_response([dict(row) for row in data])


def id_convert(from_id, to_type):
	"""Convert between different types of gene symbol (identifier) through database querying.
		Joining is used whenever possible, instead of this function.
	"""
	from_type = identifier_type(from_id)
	if from_type == to_type:
		return from_id
	sql_statement = "SELECT id2{1}.accession FROM id2{1} INNER JOIN id2{0} ON id2{1}._id=id2{0}._id WHERE id2{0}.accession=?".format(
		from_type, to_type)
	params = [from_id]
	result = db.execute(sql_statement, params).fetchall()
	result = [row[0] for row in result]
	if isinstance(result, list) and len(result) == 1:
		result = result[0]
	return result or None

def id2APA(ID, dataset):
	from_type = identifier_type(ID)
	to_type = "ensembl"
	if from_type == to_type:
		sql_statement = "SELECT d.id FROM {}_event AS d WHERE d.gene_id=?".format(dataset)
	else:
		sql_statement = "SELECT d.id FROM {2}_event AS d INNER JOIN id2{1} ON d.gene_id=id2{1}.accession INNER JOIN id2{0} ON id2{1}._id=id2{0}._id WHERE id2{0}.accession=?".format(
			from_type, to_type, dataset)
	params = [ID]
	APA_IDs = db.execute(sql_statement, params).fetchall()
	APA_IDs = [row[0] for row in APA_IDs]
	return APA_IDs

def add_trendline(fig, df, xlabel="PDUI", ylabel="TPM", init_visible=True, limited_as_ratio=(True, False)):
	"""Helper function to add a trend line to the current figure.
		Return whether fitting succeed.
	"""
	mask_labels = [label for label, limited in zip([xlabel, ylabel], limited_as_ratio) if limited]
	mask_expr = " & ".join(["0.00001 < {} < 0.99999".format(label) for label in mask_labels])	# Mask out 0s and 1s
	if mask_expr:
		df_masked = df.query(mask_expr)
	else:
		df_masked = df
	if len(df_masked) < 3:
		return False
	linear_fit_result = sm.OLS(df_masked[ylabel], sm.add_constant(df_masked[xlabel]), missing="raise").fit()
	plot_line_x = np.linspace(df[xlabel].min(), df[xlabel].max(), 21)
	plot_line_result = linear_fit_result.get_prediction(sm.add_constant(plot_line_x))
	fig.add_trace(go.Scatter(
		name="Linear regression:<br>y={:.4f}x{:+.4f}<br>r={:.4f}<br>p={:.4f}".format(
			linear_fit_result.params[xlabel],
			linear_fit_result.params["const"],
			(operator.pos if linear_fit_result.params[xlabel] > 0 else operator.neg)(math.sqrt(linear_fit_result.rsquared)), 
			linear_fit_result.pvalues[xlabel]
		),
		x=plot_line_x,
		y=plot_line_result.predicted_mean,
		mode="lines",
		hoverinfo="skip", hovertemplate=None,	# Disable line hover tip
		visible=init_visible
	))
	# Confidence interval region plot
	CI = plot_line_result.conf_int()
	interval_l = CI[:,0]
	interval_u = CI[:,1]
	fig.add_trace(go.Scatter(
		name="95% CI",
		x=np.concatenate((plot_line_x, np.flip(plot_line_x))),	# To form a closed shape
		y=np.concatenate((interval_l, np.flip(interval_u))),
		mode="none",
		fill="toself",
		fillcolor="rgba(0,0,0,0.1)",
		hoverinfo="skip", hovertemplate=None,
		visible=init_visible
	))
	return True

@app.route("/<any(QAPA,DaPars):dataset>/analysis_plot", methods=["GET"])
def analysis_plot(dataset):
	"""(This URL is invisible to user)
		Make expression-PDUI correlation plot as JSON and return it to Plotly.js
		Parameters (plotted gene) are encoded in URL as key-value pairs.
	"""
	dataset_name = dataset
	dataset = APA_datasets[dataset]

	args = request.args
	APA_ID = args["ID"]	# When error, a HTTP 400 is automatically returned
	gene_ID = args["geneID"]
	try:
		PDUI = dataset.PDUI.loc[APA_ID,:]
		TPM = dataset.TPM.loc[gene_ID,:]
	except KeyError as e:
		abort(400, description="Can't find gene ID in the dataset: \"{}\"".format(e.args[0]))

	case_type = APA_datasets["case_type"]

	df = DataFrame({"case_type": case_type.to_series().astype("category"), "PDUI":PDUI.to_series(), "TPM":TPM.to_series()})	# "case_ID":PDUI.coords["sample"]
	# pd.to_pickle(df, "test.pickle")
	df.dropna(inplace=True)
	df.reset_index(inplace=True)	# Convert index (sample) to a column
	fig = px.scatter(df, x="PDUI", y="TPM", color="case_type", labels={"case_type":"Group", "PDUI":("PAU" if dataset_name == "QAPA" else "PDUI"), "TPM":"log(TPM+1)"},
		render_mode="svg", hover_name="sample")

	add_trendline(fig, df)

	if dataset_name == "DaPars":
		# Form index for each trace (including the last two i.e. total trend line)
		n_groups = len(fig.data)-2
		indexes = list(range(1, n_groups+1))+[0, 0]
		# Make trend line for each sample group
		groups = df.groupby("case_type")
		for i, trace in enumerate(fig.data[:-2], start=1):
			name = trace.name
			succeed = add_trendline(fig, groups.get_group(name), init_visible=False)
			if succeed:
				indexes.extend([i, i])
		indexes = np.array(indexes)
		# Samples dropdown
		btns = [dict(method="restyle", label="all", args=[{"visible":[True]*(n_groups+2)+[False]*(len(indexes)-n_groups-2)}])]
		btns.extend([dict(method="restyle", label=trace.name, args=[{"visible":(indexes==i).tolist()}]) for i, trace in enumerate(fig.data[:n_groups], start=1)])
		fig.update_layout(updatemenus=[dict(
			type = "dropdown",
			x=0.07,
			xanchor="left",
			yanchor="bottom",
			showactive=False,
			pad=dict(b=5),
			buttons=btns
		)])
		fig.update_layout(annotations=[dict(
			text="Group:", showarrow=False, x=0, y=1.06, xref="paper", yref="paper", align="left"
		)])

	# fig.update_layout(legend=dict(yanchor="top", y=0.99, xanchor="right", x=0.99, bgcolor="rgba(255,255,255,0)"))
	return make_json_response(fig.to_json())

def star_symbol(p):
	"""Significance symbol under violin plot"""
	if p > 0.05:
		return "ns"
	else:
		l = -math.log10(p)
		nstar = math.floor(l)
		if nstar > 4:
			nstar = 4
		return "*"*nstar


@app.route("/<any(QAPA,DaPars):dataset>/gene_correlation_plot", methods=["POST"])
def correlation_plot(dataset):
	"""(This URL is invisible to user)
		Make expression-expression (and also PDUI-PDUI) correlation plot
		as JSON and return it to Plotly.js
		Parameters are passed in as JSON by POST.
	"""
	# Example: {"gene1":"MKI67","gene2":"TP53","selected_indexes":[1,3,4,7,8,9,...]}
	sample_names = db.execute("SELECT sample_id FROM {}_case_info ORDER BY rowid".format(dataset)).fetchall()
	sample_names = [name[0] for name in sample_names]
	data = request.get_json()	# When error, a HTTP 400 is automatically returned
	try:
		ID1 = data["gene1"]
		ID2 = data["gene2"]
		selected_indexes = data["selected_indexes"]
		selected_indexes = np.array(selected_indexes, dtype=np.int32)
	except (KeyError, TypeError, ValueError):
		abort(400, description="Submitted data structure can't be recognized")
	request_assert(isinstance(ID1, str))
	request_assert(isinstance(ID2, str))
	request_assert(selected_indexes.size >= 1)
	request_assert(selected_indexes.min() >= 1 and selected_indexes.max() <= len(sample_names))

	# Beacuse use_expression_dataset is changed to "QAPA", PDUI correlation plots can't be done anymore, and are commented out from here
	ens1 = id_convert(ID1, "ensembl")
	request_assert(isinstance(ens1, str), description="None or multiple genes are found with the specified ID")	# Single item instead of a list
	# id1_type = identifier_type(ID1)
	# convert_template = "SELECT DaPars_event.id FROM DaPars_event INNER JOIN id2refseq ON id2refseq.accession=DaPars_event.id INNER JOIN id2{0} ON id2{0}._id=id2refseq._id WHERE id2{0}.accession=?"
	# if id1_type == "refseq":
	# 	refseq1 = ID1
	# else:
	# 	# A gene corresponds to multiple transcripts (RefSeq ID), so can't use id_convert directly
	# 	refseq1 = db.execute(convert_template.format(id1_type), [ID1]).fetchall()
	# 	request_assert(len(refseq1) == 1, description="None or multiple genes are found with the specified ID")
	# 	refseq1 = refseq1[0][0]
	
	ens2 = id_convert(ID2, "ensembl")
	request_assert(isinstance(ens2, str), description="None or multiple genes are found with the specified ID")
	# id2_type = identifier_type(ID2)
	# if id2_type == "refseq":
	# 	refseq2 = ID2
	# else:
	# 	refseq2 = db.execute(convert_template.format(id2_type), [ID2]).fetchall()
	# 	request_assert(len(refseq2) == 1, description="None or multiple genes are found with the specified ID")
	# 	refseq2 = refseq2[0][0]
	# app.logger.debug("Converted IDs for correlation plot: {} {} {} {}".format(ens1, refseq1, ens2, refseq2))

	try:
		expression1 = APA_datasets[dataset].TPM.loc[ens1, sample_names]
		expression2 = APA_datasets[dataset].TPM.loc[ens2, sample_names]
		# pdui1 = APA_datasets[config["use_expression_dataset"]].PDUI.loc[refseq1, sample_names]
		# pdui2 = APA_datasets[config["use_expression_dataset"]].PDUI.loc[refseq2, sample_names]
	except KeyError as e:
		abort(400, description="Can't find gene ID in the dataset: \"{}\"".format(e.args[0]))

	case_type = APA_datasets["case_type"]
	selected_samples = np.array(sample_names)[selected_indexes-1]	# Note that selected_indexes is 1-based

	# fig = make_subplots(rows=1, cols=2, shared_xaxes=False, shared_yaxes=False, subplot_titles=["log(TPM+1)","PDUI"])

	df = DataFrame({"case_type": case_type.to_series().astype("category"), "TPM1":expression1.to_series(), "TPM2":expression2.to_series()})
	df = df.loc[selected_samples]
	df.dropna(inplace=True)
	df.reset_index(inplace=True)
	fig1 = px.scatter(df, x="TPM1", y="TPM2", color="case_type", labels={"case_type":"Group", "TPM1":"log(TPM+1) of "+ID1, "TPM2":"log(TPM+1) of "+ID2},
		render_mode="svg", hover_name="sample")
	add_trendline(fig1, df, xlabel="TPM1", ylabel="TPM2", limited_as_ratio=(False, False))
	# fig1.update_traces(legendgroup="1")
	fig1.update_yaxes(scaleanchor="x", scaleratio=1)
	# fig.add_traces(fig1.data, rows=1, cols=1)

	# df = DataFrame({"case_type": case_type.to_series().astype("category"), "PDUI1":pdui1.to_series(), "PDUI2":pdui2.to_series()})
	# df = df.loc[selected_samples]
	# df.dropna(inplace=True)
	# df.reset_index(inplace=True)
	# fig2 = px.scatter(df, x="PDUI1", y="PDUI2", color="case_type", labels={"case_type":"Group", "PDUI1":ID1+" PDUI", "PDUI2":ID2+" PDUI"},
	# 	render_mode="svg", hover_name="sample")
	# fig2.update_traces(showlegend=False)
	# add_trendline(fig2, df, xlabel="PDUI1", ylabel="PDUI2", limited_as_ratio=(True, True))
	# fig2.update_traces(legendgroup="2")
	# # fig2.update_yaxes(scaleanchor="x", scaleratio=1)
	# fig.add_traces(fig2.data, rows=1, cols=2)

	# FIXME: Any name conflict of two plots?
	# BUG: update_yaxes() will lead to obviously wrong result
	# fig.update_layout(legend=dict(yanchor="top", y=0.99, xanchor="right", x=0.99, bgcolor="rgba(255,255,255,0)"))
	# return make_json_response(fig.to_json())

	return make_json_response(fig1.to_json())

@app.route("/<any(QAPA,DaPars):dataset>/de_analysis_plot", methods=["POST"])
def DE_analysis_plot(dataset):
	"""(This URL is invisible to user)
		Make violin plot as JSON and return it to Plotly.js
		Parameters are passed in as JSON by POST.
	"""
	# Example: {"gene":"MKI67","group_name_indexes":{"normal":1,"patient":2},"sample_groups":[null,null,2,1,2,null,...]}
	sample_names = db.execute("SELECT sample_id FROM {}_case_info ORDER BY rowid".format(dataset)).fetchall()
	sample_names = [name[0] for name in sample_names]
	data = request.get_json()
	# TODO: Validate request JSON
	try:
		ID = data["gene"]
		group_name_index = data["group_name_indexes"]
		sample_groups = data["sample_groups"]
	except KeyError:
		abort(400, description="Submitted data structure can't be recognized")
	request_assert(isinstance(group_name_index, dict))
	request_assert(all(isinstance(k, str) and isinstance(v, int) for k,v in group_name_index.items()))
	request_assert(len(sample_groups) == len(sample_names))
	request_assert(all(i is None or isinstance(i, int) for i in sample_groups))

	ens = id_convert(ID, "ensembl")
	request_assert(isinstance(ens, str), description="None or multiple genes are found with the specified ID")
	# APA_datasets["DaPars"].TPM is actually feature_counts

	APA_IDs = id2APA(ID, dataset)

	try:
		expression = APA_datasets[dataset].TPM.loc[ens, sample_names]
		per_site_expression = None
		if dataset == "QAPA":
			per_site_expression = APA_datasets["QAPA"].site_TPM.loc[APA_IDs, sample_names]
		per_site_PDUI = APA_datasets[dataset].PDUI.loc[APA_IDs, sample_names]
	except KeyError as e:
		abort(400, description="Can't find gene ID in the dataset: \"{}\"".format(e.args[0]))
	
	sample_groups = Series(sample_groups, index=sample_names, dtype="category", name="group")
	sample_groups = sample_groups.cat.rename_categories({v:k for k,v in group_name_index.items()})
	df = DataFrame({"group":sample_groups, "all":expression.to_series()})
	if dataset == "QAPA":
		per_site_expression = per_site_expression.transpose().to_pandas()
		per_site_expression.rename(columns=lambda site:site.split("_",1)[1], inplace=True)
		df = df.join(per_site_expression)
	per_site_PDUI = per_site_PDUI.transpose().to_pandas()
	suffix = "_PAU" if dataset == "QAPA" else "_PDUI"
	per_site_PDUI.rename(columns=lambda site:site.split("_",1)[1]+suffix, inplace=True)
	df = df.join(per_site_PDUI)
	df.index.name = "sample"
	# pd.to_pickle(df, "test.pickle")

	if dataset == "QAPA":
		n_cols = (df.shape[1]-1)//2
		width_ratio = 0.7 / (n_cols+1)
		fig = make_subplots(rows=1, cols=3, horizontal_spacing=0.02,
			column_widths=[width_ratio, 0.7-width_ratio, 0.3], y_title="log(TPM+1)")# shared_yaxes=True, x_title="APA isoform", )
	else:
		fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.02,
			column_widths=[0.7, 0.3], y_title="log(TPM+1)")# shared_yaxes=True, x_title="APA isoform", )
	groups = df.groupby("group")
	groups = [(l,g.drop(columns="group")) for l,g in groups]
	for i, (group_name, g) in enumerate(groups):
		for col_name, col in g.items():
			if col_name.endswith(suffix):
				# Ridgeline plot of sample PDUIs
				fig.add_trace(go.Violin(x=col, y=[col_name]*len(col), orientation="h", side="positive", marker_color=plotly_colors[i%len(plotly_colors)], meanline_visible=True, points=False, box_visible=False,
					legendgroup=group_name, scalegroup=group_name+"_"+col_name, offsetgroup=group_name, showlegend=False), row=1, col=3 if dataset == "QAPA" else 2)
			else:
				fig.add_trace(go.Violin(x=[col_name]*len(col), y=col, name=group_name, box_visible=True, hoveron="violins+points", marker_color=plotly_colors[i%len(plotly_colors)],
					legendgroup=group_name, scalegroup=group_name+"_"+col_name, offsetgroup=group_name, showlegend=(col_name == "all")), row=1, col=(1 if col_name == "all" else 2))

	groups = [g for l,g in groups]
	annotations = {}
	if len(data["group_name_indexes"]) <= 1:
		abort(400)
	elif len(data["group_name_indexes"]) == 2:
		# t test
		group1, group2 = groups
		for (i1,g1), (i2,g2) in zip(group1.items(), group2.items()):
			assert i1 == i2
			g1.dropna(inplace=True)
			g2.dropna(inplace=True)
			tstat, p, dof = ttest_ind(g1, g2)
			if math.isfinite(p):
				tag = "<b>{}</b><br>(p={:.4f})".format(star_symbol(p), p)
				if i1.endswith(suffix):
					n_last_subplot = 3 if dataset == "QAPA" else 2
					fig.add_annotation(y=i1, yref="y{}".format(n_last_subplot), x=0, xref="x{} domain".format(n_last_subplot),
						showarrow=False, text=tag)
				else:
					fig.add_annotation(x=i1, xref="x" if i1 == "all" else "x2", y=0, yref=("y" if i1 == "all" else "y2")+" domain",
						showarrow=False, text=tag)
	else:
		# ANOVA
		for cols in zip(*[g.items() for g in groups]):
			names = [c[0] for c in cols]
			values = [c[1] for c in cols]
			assert all(n == names[0] for n in names[1:])
			result = anova_oneway([a.dropna() for a in values], use_var="equal")
			if math.isfinite(result.pvalue):
				tag = "<b>{}</b><br>(p={:.4f})".format(star_symbol(result.pvalue), result.pvalue)
				if names[0].endswith(suffix):
					n_last_subplot = 3 if dataset == "QAPA" else 2
					fig.add_annotation(y=names[0], yref="y{}".format(n_last_subplot), x=0, xref="x{} domain".format(n_last_subplot),
						showarrow=False, text=tag)
				else:
					fig.add_annotation(x=names[0], xref="x" if names[0] == "all" else "x2", y=0, yref=("y" if names[0] == "all" else "y2")+" domain",
						showarrow=False, text=tag)

	fig.update_layout(violinmode="group")
	# fig.update_yaxes(rangemode="tozero")
	fig.update_layout(violingap=0.15, violingroupgap=0.1)
	fig.update_traces(box_width=0.12)
	# fig.update_xaxes(matches=None)
	return make_json_response(fig.to_json())


@app.route("/<any(QAPA,DaPars):dataset>/survival_plot", methods=["POST"])
def survival_plot(dataset):
	"""(This URL is invisible to user)
		Make survival plot as JSON and return it to Plotly.js
		Parameters are passed in as JSON by POST, and the data format is similar to
		de_analysis_plot (when specifying groups) or correlation_plot (when specifying a gene).
	"""

	data = request.get_json()
	# The same as in de_analysis_plot
	try:
		survival_variable = data["survival_variable"]
		if survival_variable == "group":
			group_name_index = data["group_name_indexes"]
			sample_groups = data["sample_groups"]
			request_assert(isinstance(group_name_index, dict))
			request_assert(all(isinstance(k, str) and isinstance(v, int) for k,v in group_name_index.items()))
			request_assert(all(i is None or isinstance(i, int) for i in sample_groups))
		elif survival_variable == "gene":
			selected_indexes = data["selected_indexes"]
			selected_indexes = np.array(selected_indexes, dtype=np.int32)
			ID = data["gene"]
		else:
			abort(400, description="Type of variable in survival analysis can't be recognized")
	except (KeyError, TypeError, ValueError):
		abort(400, description="Submitted data structure can't be recognized")

	df = pd.read_sql_query("SELECT vital_status, days_to_death, days_to_last_followup FROM {}_case_info ORDER BY rowid".format(dataset), db)
	df["is_dead"] = df["vital_status"] == "Dead"
	df["alive_duration"] = df["days_to_death"].where(df["is_dead"], other=df["days_to_last_followup"])
	df = df[["is_dead", "alive_duration"]]

	
	if survival_variable == "gene":
		ID = id2APA(ID, dataset)
		if len(ID) != 1:
			abort(400, description="None or multiple genes are found")
		ID = ID[0]
		pdui = APA_datasets[dataset].PDUI.loc[ID,:].to_series().reset_index(drop=True)
		threshold = np.nanmedian(pdui)
		if threshold == 0:
			if not np.any(pdui>threshold):
				abort(400, description="Too many samples have PDUI=0.")
			else:
				sample_groups = pdui>threshold
		else:
			if threshold == 1 and not np.any(pdui<threshold):
				abort(400, description="Too many samples have PDUI=1.")
			else:
				sample_groups = pdui>=threshold
		sample_groups = Series(sample_groups, dtype="category")
		sample_groups = sample_groups.cat.rename_categories({False:"low PDUI",True:"high PDUI"})
		df = df.assign(group=sample_groups)
		df = df.iloc[selected_indexes-1]	# Note that selected_indexes is 1-based
	else:
		sample_groups = Series(sample_groups, dtype="category")
		sample_groups = sample_groups.cat.rename_categories({v:k for k,v in group_name_index.items()})
		df = df.assign(group=sample_groups)
	df.dropna(inplace=True)
	chisq, p = sm.duration.survdiff(df["alive_duration"], df["is_dead"], df["group"])
	groups = df.groupby("group")
	survfuncs = [sm.SurvfuncRight(d["alive_duration"], d["is_dead"], title=str(k)) for k,d in groups]
	anno = ["MST of {}={}".format(s.title, s.quantile(0.5)) for s in survfuncs]	# Median survival time
	anno.append("log-rank p={:.4f}".format(p))

	fig = go.Figure()
	plot_survfunc(survfuncs, fig)
	anno = "<br>".join(anno)
	fig.add_annotation(text=anno, x=0.1, y=0.95, xref="x domain", yref="y domain", showarrow=False)
	fig.update_xaxes(title_text='Time')
	fig.update_yaxes(title_text='Proportion Survival')
	return make_json_response(fig.to_json())

def plot_survfunc(survfuncs, fig):
	"""From statsmodels.duration.survfunc.plot_survfunc(). With modification to use Plotly.
	"""

	if isinstance(survfuncs, sm.SurvfuncRight):
		survfuncs = [survfuncs]
	else:
		assert all(isinstance(s, sm.SurvfuncRight) for s in survfuncs)

	for gx, sf in enumerate(survfuncs):
		# The estimated survival function does not include a point at
		# time 0, include it here for plotting.
		surv_times = np.concatenate(([0], sf.surv_times))
		surv_prob = np.concatenate(([1], sf.surv_prob))

		# If the final times are censoring times they are not included
		# in the survival function so we add them here
		mxt = max(sf.time)
		if mxt > surv_times[-1]:
			surv_times = np.concatenate((surv_times, [mxt]))
			surv_prob = np.concatenate((surv_prob, [surv_prob[-1]]))

		label = getattr(sf, "title", "Group %d" % (gx + 1))

		fig.add_trace(go.Scatter(x=surv_times, y=surv_prob, name=label,
			mode="lines", line_shape="hv", marker_color=plotly_colors[gx%len(plotly_colors)]))

		# Plot the censored points.
		ii = np.flatnonzero(np.logical_not(sf.status))
		ti = np.unique(sf.time[ii])
		jj = np.searchsorted(surv_times, ti) - 1
		sp = surv_prob[jj]
		fig.add_trace(go.Scatter(x=ti, y=sp, name="censored", mode="markers",
			marker_symbol="cross-thin-open", marker_color=plotly_colors[gx%len(plotly_colors)]))

	fig.update_layout(yaxis_range=[0, 1.01])


@app.route("/<any(QAPA,DaPars):dataset>/expression_summary_plot", methods=["GET"])
def expression_summary_plot(dataset):
	"""(This URL is invisible to user)
		Make expression box + points plot
		as JSON and return it to Plotly.js
		Parameters (gene symbol) are encoded in URL as key-value pairs.
	"""
	dataset = APA_datasets[dataset]

	args = request.args
	site_ID = args["siteID"]	# When error, a HTTP 400 is automatically returned
	try:
		TPM = dataset.site_TPM.loc[site_ID,:]
	except KeyError as e:
		abort(400, description="Can't find gene ID in the dataset: \"{}\"".format(e.args[0]))
	case_type = APA_datasets["case_type"]

	df = DataFrame({"case_type": case_type.to_series().astype("category"), "TPM":TPM.to_series()})
	df.dropna(inplace=True)
	fig = px.box(df, x="case_type", y="TPM", labels={"case_type":"Group", "TPM":"log(TPM+1)"}, points="all")
	fig.update_traces(pointpos=0, jitter=0.6, marker_opacity=0.5)
	# fig.update_layout(boxgap=0)
	return make_json_response(fig.to_json())


@app.route("/gene_detail", methods=["GET"])
def gene_detail():
	"""(This URL is invisible to user)
		Query and return gene information as JSON to the popover.
		Query parameters (gene symbol) are encoded in URL as key-value pairs.
	"""
	ID = request.args["ID"]	# When error, a HTTP 400 is automatically returned
	ID_type = identifier_type(ID)
	request_assert(ID_type != "symbol")
	# Currently I don't know how to retrieve these columns in a single query. (INNER JOIN will produce extra rows because of Cartesian product.)
	data = db.execute("SELECT d1._id,d2.accession,d2.full_name,d2.gene_type FROM id2symbol AS d2 INNER JOIN id2{} AS d1 ON d1._id=d2._id WHERE d1.accession=?".format(ID_type), [ID]).fetchall()
	request_assert(len(data) == 1, description="None or multiple genes are found with the specified ID")
	data = dict(data[0])
	entrez = db.execute("SELECT accession AS entrez FROM id2entrez WHERE _id=?", [data["_id"]]).fetchall()
	assert len(entrez) == 1
	entrez = entrez[0]
	alias = db.execute("SELECT group_concat(accession) AS alias FROM id2alias WHERE _id=? GROUP BY _id", [data["_id"]]).fetchall()
	assert len(alias) == 1
	alias = alias[0]
	band = db.execute("SELECT group_concat(band) AS band FROM id2band WHERE _id=? GROUP BY _id", [data["_id"]]).fetchall()
	assert len(band) == 1
	band = band[0]
	location = db.execute("SELECT group_concat(printf('%s:%d-%d',chrom,gene_start,gene_end),'<br>') AS location FROM id2location WHERE _id=?", [data["_id"]]).fetchall()
	location = location[0]
	data.update(entrez)
	data.update(alias)
	data.update(band)
	data.update(location)
	return make_json_response(data)


@app.errorhandler(BadRequest)
def handle_BadRequest(e):
	"""Return JSON instead of HTML for HTTP 400 (Bad Request).
		Now all pages which need args return JSON.
	"""
	app.logger.debug("400 error description: "+e.description)
	return jsonify(code=e.code, name=e.name, description=e.description), 400
