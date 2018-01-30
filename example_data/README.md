This is an example dataset for testing. The data are a replication of the experiment described in:

Macklaim, Jean M., et al. "Comparative meta-RNA-seq of the vaginal microbiota and differential expression by Lactobacillus iners in health and dysbiosis." Microbiome 1.1 (2013): 12.

The primary comparison of these data are BV status (by Nugent scoring).
- BV samples: "012B","013A","014B","017B","018B"
- non-BV (healthy) samples: "010B", "001B", "009B"

The counts table contains 8 samples (columns) total, with the number of reads mapped per sample to each reference subsys4 (rows)

ALL samples have been mapped to Jean's custom SEED subsys4 reference set (curated from http://www.theseed.org) to generate the counts table. The samples are columns, with the subsys4 features as rows. Note that after subsys4, the functional assignments are hierachrical and overlapping, whereas each subsys4 bin is unique.

For reading into R:
`d<-read.table("counts.txt", header=T, sep="\t", quote="", check.names=F, row.names=1)`
