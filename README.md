# SRIQ
This is a implementation of Systematic Random forest Integrative Qualitative threshold (SRIQ) clustering

To apply this pipe-line on your gene expression data, clone this repository to desired folder 

>The format your expression data should be in is:
>First row: Gene \t Observation1 \t Observation2 \t ... \t ObservationN
>Second row: gene1 \t val1 \t val2 \t ... \t valN
>Third row: gene2 \t val1 \t val2 \t ... \t valN
>...
>N:th row: gene(N-1) \t val1 \t val2 \t ... \t valN

Pipeline assumes that the data is fpkm normalized.
