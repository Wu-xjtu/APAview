# APAview

## Introduction

**APAview** is a web-based platform for **alternative polyadenylation** (APA) features exploration on APA quantification data from QAPA/DaPars in hematological tumors. Taking APA, gene/transcript expression, and clinical data as input, APAview could:

**(1)** provide correlation analysis between APA usage index and gene/transcript expression;

**(2)** provide differential analysis of APA usage index or other features among groups; 

**(3)** identify genes with shortened/lengthened 3’UTR between groups; 

**(4)** provide survival analysis based on APA usage index;

**(5)** annotate genes and APA sites using databases such as NCBI and polyA DB;

**(6)**  visualize gene structures, APA sites, miRNA, RBP, and APA motifs in IGV plots.

## Installation and Usage

1. Install Python 3.9 on Windows/Linux.

2. Create  an environment for APAview running.

   Windows:

   ```
   python -m venv APAview_env
   .\APAview_env\Scripts\Activate.ps1 (using python)
   .\APAview_env\Scripts\activate.bat (using cmd)
   ```

   Linux:

   ```
   python3 -m venv APAview_env
   . ./APAview_env/bin/activate
   ```

3. Install dependency packages.

   Windows:

   ```
   python -m pip install --upgrade setuptools pip
   python -m pip install wheel
   python -m pip install flask statsmodels xarray h5netcdf plotly
   #########################
   Package         Version
   --------------- --------
   Flask           2.0.3
   h5netcdf        0.14.1
   pandas          1.4.1
   plotly          5.6.0
   statsmodels     0.13.2
   xarray          2022.3.0
   ##########################
   ```

   Linux:

   ```
   python3 -m pip install --upgrade setuptools pip
   python3 -m pip install wheel
   python3 -m pip install flask statsmodels xarray h5netcdf plotly
   #########################
   Package         Version
   --------------- --------
   Flask           2.0.3
   h5netcdf        0.14.1
   pandas          1.4.1
   plotly          5.6.0
   statsmodels     0.13.2
   xarray          2022.3.0
   ##########################
   ```

4. Download the latest development version of APAview and unzip the file to `/pathway/APAview`.

5. Run APAview in `APAview_env`.

   ```
   cd /pathway/APAview
   flask run
   #################
    * Environment: production
      WARNING: This is a development server. Do not use it in a production deployment.
      Use a production WSGI server instead.
    * Debug mode: off
    * Running on http://127.0.0.1:5000/ (Press CTRL+C to quit)
    ################
   ```

6. Open the link (http://127.0.0.1:5000/) of APAview on the browser (e.g. Firefox).

## Input data format

1. APA quantification data

   1.1 QAPA output (*https://www.github.com/morrislab/qapa*)

   | Column         | Description                                                  |
   | :------------- | ------------------------------------------------------------ |
   | APA_ID         | unique identifier in the format "<Ensembl Gene ID> _ <number> _ <(P\|D\|S)>", P = proximal, D = distal, and S = single |
   | Transcript     | one or more Ensembl Transcript IDs                           |
   | Gene           | Ensembl Gene ID                                              |
   | Gene_Name      | gene symbol                                                  |
   | Chr            | chromosome                                                   |
   | LastExon.Start | start coordinate of last exon                                |
   | LastExon.End   | end coordinate of last exon                                  |
   | Strand         | "+" or "-"                                                   |
   | UTR3.Start     | start coordinate of 3'UTR                                    |
   | UTR3.End       | end coordinate of 3'UTR                                      |
   | Length         | length of the 3′UTR                                          |
   | Num_Events     | number of PAS per gene                                       |
   | <sample1>.PAU  | PolyA site Usage (PAU) estimate for sample1                  |
   | <sample1>.TPM  | TPM estimate for sample1                                     |

   *example:*

   ![image-20220328195442643](C:\Users\lenovo\AppData\Roaming\Typora\typora-user-images\image-20220328195442643.png)

   1.2 DaPars output (*http://bioinfo.szbl.ac.cn/DaPars2/DaPars2.html*)

   | Column                 | Description                                                  |
   | ---------------------- | ------------------------------------------------------------ |
   | Gene                   | transcripts information in "RefSeq_transcript_ID\| gene_symbol \| chromosome \| strand" |
   | Predicted_Proximal_APA | predicted proximal APA site by DaPars                        |
   | Loci                   | the 3'UTR region of the transcript                           |
   | Other columns          | the Percentage of Distal APA site Usage Index (PDUI) values of each sample calculated by DaPars |

   *example:*

   ![image-20220328194812947](C:\Users\lenovo\AppData\Roaming\Typora\typora-user-images\image-20220328194812947.png)

2. gene/transcription data

   Rows: genes/transcripts using gene symbols/IDs as rows' names.

   Columns: expression values quantified in *TPM* format, each column represents a sample.

   *example:*

   <img src="C:\Users\lenovo\AppData\Roaming\Typora\typora-user-images\image-20220328195707846.png" alt="image-20220328195707846" style="zoom:80%;" />

3. clinical data

   Rows: each row represents a sample.

   Columns: information like age, gender, disease types, survival data, and so on. The first column is *sample IDs* consistent with that in APA and expression data. Survival information should contain three columns including *vital_status*, *days_to_last_followup*, and *days_to_death*.

   *example:*

   <img src="C:\Users\lenovo\AppData\Roaming\Typora\typora-user-images\image-20220328194618437.png" alt="image-20220328194618437" style="zoom: 80%;" />

## Data import

1.  Put all imported files into a directory (e.g. `source_dir`).
2.  Run `import_data.sh source_dir/` in APAview_env.
3.  Generate `APA.sqlite` and `data.nc` in the current directory. The two files store all input information. 

## Citation

*APAview: a web-based platform for alternative polyadenylation analyses in hematological tumors*

## Contact

Emails for any comments, suggestions, questions, and bug reports of APAview are welcomed to send to Xiaoming Wu([wxm@mail.xjtu.edu.cn](mailto:wxm@mail.xjtu.edu.cn)).