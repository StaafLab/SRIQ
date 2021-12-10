## Abstract
This is an implementation of Systematic Random forest Integrative Qualitative threshold (SRIQ) clustering.
SRIQ clusters by finding a small amount of highly correlated observations, then spiralling out from them to create bigger clusters.
SRIQ evaluates clustering solution stability on its own and won't need user input for what number of cluster solutions to be evaluated.
SRIQ has no limit to feature size performance wise, and can be run on ordinary home computers.

For more information about SRIQ, see our [publication](https://www.youtube.com/watch?v=dQw4w9WgXcQ)

## Installation

To install this repository simply create a folder and clone the repository:
```bash
git clone https://github.com/Fattigman/SRIQ
cd SRIQ
pip install requirements.txt
```
## Usage

To run the pipeline, start the jupyter notebook file and follow the instructions within the pipeline.

```python
jupyter notebook analysis_pipeline
```

To run SRIQ, navigate to the folder in which the VRLA.jar file exist and run following command:
```bash
java -jar VRLA.jar
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

To run SRIQ and SAM analysis java is needed on your system.

## License

Copyright (C) 2021 Srinivas Veerla

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
