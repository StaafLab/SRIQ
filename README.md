## Abstract
This is an implementation of Systematic Random forest Integration to Qualitative threshold (SRIQ) clustering.
SRIQ clusters by finding a core clusters of highly correlated observations, then spiralling out from them to create bigger clusters.
SRIQ evaluates clustering solution stability on its own and won't need user input for what number of cluster solutions to be evaluated.
SRIQ has no limit to feature size performance wise, and can be run on ordinary home computers.

For more information about SRIQ, see our [publication](https://www.youtube.com/watch?v=dQw4w9WgXcQ)


## To Run SRIQ (JAVA)
<b>Step 1:</b> Data Pre-process,navigate to the folder in which the SRIQPreprocess.jar file exist and run following command:
```bash
java -jar SRIQPreprocess.jar path-to/datafile.txt
```
<b>output:</b> a file with name ending with "_mc_log2_nz" will be in the same path folder as input file folder.<br><br>
<b>Step 2:</b> Edit the test.properties file<br>
#testing
#Fri Mar 15 14:39:14 CET 2019
studyName=LU_LUAD_SRIQ
studyPath=/home/Researcher/SRIQ/software/data/
inFileName=GeneData_AC_fpkm_mc_log2_nz
outPath=/home/Researcher/SRIQ/software/output/
distCutOff=0.8, 0.79, 0.78, 0.77, 0.76, 0.75, 0.74, 0.73, 0.72, 0.71, 0.7, 0.69, 0.68, 0.67, 0.66, 0.65, 0.64, 0.63, 0.62, 0.61, 0.6, 0.59, 0.58, 0.57, 0.56, 0.55, 0.54, 0.53, 0.52
permutations=10000
iterations=10
spiral=TRUE
minClusterSize=0
minBagSize=1200
method=PEARSON




## Installation

To install this repository simply create a folder and clone the repository:
```bash
git clone https://github.com/StaafLab/SRIQ
cd SRIQ
pip install -r requirements.txt
```

## Usage

To run the pipeline, start the jupyter notebook file and follow the instructions within the pipeline.
It's recommended to use virtual environments when running the pipeline in order to avoid conflicting packages.

```python
jupyter notebook analysis_pipeline
```

To run SRIQ, navigate to the folder in which the SRIQ.jar file exist and run following command:
```bash
java -jar SRIQ.jar path-to/test.properties
```
Alternative to run SAMDEG, navigate to the folder in which the SAMDEG.jar file exist and run following command:
```bash
java -jar <path-to/SAMDEG.jar> <path-to/test.properties> <spiral (true or false)> <diameter> <no. of clusters> <q-value> <fold-change> <log2_transformed_gex_file>
e.g., java -jar SAMDEG.jar "F:/test/LUAD_test/test.properties" false 0.63 6 0 2 "F:/test/LUAD_test/newFiltered_35k.txt"
```
## Features

The project provide following features:

* SRIQ-clustering
* Pre-clustering data normalization
* Silhoutte-plot analysis
* UMAP
* Differential gene expression with SAM
* Gene enrichment analysis with EnrichR against customizable databases
* Visualization of single or multiple genes across clusters

## Data format requirements

Your data, to be clustered, should be a tab separated .txt file and look like this:

For SRIQ to accept the data to be clustered, the file has to be in following format:
| Gene               |   Sample1          |   Sample2          |   ...              |   SampleN          |
|:-------------------|-------------------:|-------------------:|-------------------:|-------------------:|
| Genename 1         |            val1    |         val2       |         ...        |         valN       |
| ...                |            ...     |         ...        |         ...        |         ...        |
| Genename M         |            val1    |         val2       |         ...        |         valN       |

For clustering, expression data should be:
1. Off-set by 0.1 OR set all values < 1 set to 1
2. Log2transformed
3. Median-centered

For SAM-analysis:
1. Off-set by 0.1 OR all values < 1 set to 1
2. Log2transformed

Before running SRIQ, test.properties file need to be correctly configured.
The following lines needs to be correct otherwise SRIQ won't start.
```
studyName: Desired output name for project
studyPath: Folderpath to expression data
inFileName: expression file. exclude '.txt' from file name
outPath: Folderpath for SRIQ output
```

The enrichR module assumes that gene names are in the form of gene symbols. I have implemented [mygene](https://mygene.info) api, set the variable 'scopes' as the format of your gene names.

## Package requirements
For python packages run following command:
```bash
pip install -r requirements.txt
```

To run SRIQ and SAM analysis, java is needed on your system.

## License

Copyright (C) 2021 Srinivas Veerla

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
