# SRIQ
Author = Jacob Karlström, Srinivas Veerla
Maintainer Jacob

## Abstract
This is a implementation of Systematic Random forest Integrative Qualitative threshold (SRIQ) clustering.
SRIQ clusters by finding a small amount of highly correlated observations, then spiralling out from them to create bigger clusters.
SRIQ evaluates clustering solution stability on its own and won't need user input for what number of cluster solutions to be evaluated.
SRIQ has no limit to feature size performance wise, and can be run on ordinary home computers.
For more information about SRIQ, see our [publication](https://www.youtube.com/watch?v=dQw4w9WgXcQ)

To apply this pipe-line on your gene expression data, clone this repository to desired folder 

The format your expression data should be in is:
```
First row: Gene \t Observation1 \t Observation2 \t ... \t ObservationN

Second row: gene1 \t val1 \t val2 \t ... \t valN

Third row: gene2 \t val1 \t val2 \t ... \t valN

...

N:th row: gene(N-1) \t val1 \t val2 \t ... \t valN

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
