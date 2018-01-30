## Diamond pipeline to map reads to SEED database

##### Jul 19 2017
By: Jean Macklaim

---
DIAMOND (https://github.com/bbuchfink/diamond) is a super fast sequencing aligner for protein or DNA translation (nt->protein) searches. The database was custom curated from SEED (http://www.theseed.org/wiki/Main_Page) by Jean to contain non-redundant protein sequences with assigned SEED usbsystem functions.

The SEED database:  
Briefly, each `fig.peg` sequence comes from an individual genome and has a subsys4 (enzymatic) functional assignment. There are many `fig.pegs` in a subsys4 category. From there, subsys4 are hierarchically organized into broader functional groups (subsys3, subssy2, subsys1) which are non-unique to a subsys4 (i.e. a subsys4 can belong to multiple subsys3).

As Jean for more information....

---
### Required to run
```
diamond_to_seed.sh
merge_counts.pl
blast_to_counts.pl
```

On Agrajag:
DIAMOND in installed in `/Volumes/data/bin/diamond`
Scripts installed in `/Volumes/bin/diamond-program/` (check with github version for updates)
The SEED database (pre-formatted for DIAMOND) is `/Volumes/data/SEED_database/subsys4.dmnd`

### Running a DIAMOND search

There is a shell script `diamond_to_seed.sh`. To run, do the following
- You should have a data directory (usually called `data`) containing the **demultiplexed** read files in `.fastq.gz` format/extension. The file names should look like: `F8G-2_S43_R1_001.fastq.gz` where everything before the file extension `.fastq.gz` will be taken as the sample name for downstream output
- Make a working directory e.g. `mkdir map_seed`. Copy `diamond_to_seed.sh` and the `bin` directory containing two perl scripts `blast_to_counts.pl` and `merge_counts.pl` to this directory
- Before running `diamond_to_seed.sh` define the paths at the top of script

>>Example directory structure:
>>- ../project_name/data/reads.fastq.gz
>>- ../project_name/map_seed/diamond_to_seed.sh
>>- ../project_name/map_seed/bin/blast_to_counts.pl

Running:
`nohup ./diamond_to_seed.sh &`
*Note: Since this program takes a while to run, using `nohup` will push it to the background and push any output to terminal to `nohup.out`*

### Output:
A `diamond_output` directory will be created in the directory you run the script with the .daa files, the converted .m8 files, and the total counts table `all_counts.txt` will be in your working directory

### Troubleshooting
You may need to delete the output directory and its contents to run the script again
`rm -R diamond_output`
