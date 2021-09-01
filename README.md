## Abstract
This is an implementation of Systematic Random forest Integrative Qualitative threshold (SRIQ) clustering.
SRIQ clusters by finding a small amount of highly correlated observations, then spiralling out from them to create bigger clusters.
SRIQ evaluates clustering solution stability on its own and won't need user input for what number of cluster solutions to be evaluated.
SRIQ has no limit to feature size performance wise, and can be run on ordinary home computers.

For more information about SRIQ, see our [publication](https://www.youtube.com/watch?v=dQw4w9WgXcQ)

## Installation

To install this repository simply create a folder and clone the repository:
```
mkdir SRIQ
cd SRIQ
git clone https://github.com/Fattigman/SRIQ
```

## Data requirements

The format your expression data should be in is:
```
'|    | Gene               |   TCGA-05-4384-01A |   TCGA-05-4390-01A |   TCGA-05-4396-01A |   TCGA-05-4405-01A |   TCGA-05-4410-01A |\n|---:|:-------------------|-------------------:|-------------------:|-------------------:|-------------------:|-------------------:|\n|  0 | ENSG00000242268.2  |            0.12364 |         0          |         0.148773   |         0          |            0       |\n|  1 | ENSG00000270112.3  |            0       |         0.00586692 |         0.00688089 |         0.00639122 |            0       |\n|  2 | ENSG00000167578.15 |            4.04325 |         2.0436     |         2.02117    |         3.11505    |            5.08969 |\n|  3 | ENSG00000273842.1  |            0       |         0          |         0          |         0          |            0       |\n|  4 | ENSG00000078237.5  |            4.91725 |         4.8264     |         3.28451    |         4.13261    |            4.65857 |'

```
Pipeline assumes that the data is fpkm normalized beforehand.

Before running SRIQ correct settings has to be set for SRIQ in software/VRLA/resources/test.properties.

```
studyName: Desired output name for project
studyPath: Folderpath to expression data
inFileName: txt expression file. exclude .txt values tab separated
outPath: Folderpath for SRIQ output
```
To run GO enrichment, the significantly expressed gene names need to be converted into symbol names.
To do that, run the `obj.ens2symbol(arg)`. Parameter can be found in [this documentation](https://docs.mygene.info/en/latest/doc/query_service.html)
