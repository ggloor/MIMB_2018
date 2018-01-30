Step 1) Pull genome files (see `build_database_new.sh`)

Step 2) Cluster ORFs to remove redundancy (see `cluster.sh`)

NOTES:

This workflow came from:
`/Volumes/rhamnosus/reference_genomes/may2013/cluster/uclust/`

The `.ffn` reference sequences were concatenated into one file

Because of memory limitations, I couldn't use true centroid clustering. I sorted the sequences be deviation from median length and clustered from there
