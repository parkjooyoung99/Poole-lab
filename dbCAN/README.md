# :bangbang: All work should be done in ```/workdir/XX/```:bangbang:
## 1. Load the env and Download the DB
#### 1-1. Load the env
```
source /programs/miniconda3/bin/activate run_dbcan4
export PATH=/programs/seqtk:$PATH
```
#### 1-2. Download the DB
```
test -d db || mkdir db
cd db \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/fam-substrate-mapping-08012023.tsv && mv fam-substrate-mapping-08012023.tsv fam-substrate-mapping.tsv \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/PUL.faa && makeblastdb -in PUL.faa -dbtype prot \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_12-12-2023.xlsx && mv dbCAN-PUL_12-12-2023.xlsx dbCAN-PUL.xlsx \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL.tar.gz && tar xvf dbCAN-PUL.tar.gz && rm dbCAN-PUL.tar.gz \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN_sub.hmm && hmmpress dbCAN_sub.hmm \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/CAZyDB.07262023.fa && mv CAZyDB.07262023.fa CAZyDB.fa  && diamond makedb --in CAZyDB.fa -d CAZy \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/dbCAN-HMMdb-V12.txt && mv dbCAN-HMMdb-V12.txt dbCAN.txt && hmmpress dbCAN.txt \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-1.hmm && hmmpress tf-1.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-2.hmm && hmmpress tf-2.hmm \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/stp.hmm && hmmpress stp.hmm \
    && cd ../ && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.fna \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.faa \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.gff
```

## 2. Input preparation 
:heavy_exclamation_mark: Input fasta should already be cleaned - which means all the QC are done 
:heavy_exclamation_mark: Here we use kneaddata cleaned fastq generated during the [biobakery](https://github.com/parkjooyoung99/Poole-lab/tree/b844ccc705f5eb6c7e74c349d4c632eea06a478c/Biobakery_workflow)

#### 2-1. Link fastq files
```
mkdir -p /workdir/XXX/run_dbcan/input

while read line; do echo "$line" && ln -s "/workdir/XXX/output_ours_baseline/kneaddata/main/$line" "/workdir/XXX/run_dbcan/input/$line"; done < fastqlist.txt
```


:heavy_exclamation_mark: Methods the dbCAN uses are describe [here](https://dbcan.readthedocs.io/en/latest/user_guide/run_from_raw_reads.html)
#### 2-1. 
:heavy_exclamation_mark: Uses all three methods
:small_red_triangle_down:Name should be 'Sample.R1.fastq.gz' or 'Sample.R2.fastq.gz' format
make new directory that has all input fastq.gz file
```
mkdir input_ours_all
ln -s "originalfastqpath/originalfastqfilename.gz" "/workdir/XX/input_ours_all/newfastqfilename.gz"
```
#### 2-2. Run wgmx ####
```
export BIOBAKERY_WORKFLOWS_DATABASES=/workdir/jp2626/biobakery_workflow_databases
biobakery_workflows wmgx --input input_ours_all --output output_ours_all --bypass-strain-profiling --local-jobs 8 --threads 8
```


