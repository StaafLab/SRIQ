# SRIQ
This is a implementation of Systematic Random forest Integrative Qualitative threshold (SRIQ) clustering.
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

Before running SRIQ correct settings has to be set for SRIQ in test.properties.

```
studyName: Desired output name for project
studyPath: Folderpath to expression data
inFileName: expression file. exclude .txt and separator should be \t
outPath: Folderpath for SRIQ output
```
