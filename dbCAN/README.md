# :bangbang: All work should be done in ```/workdir/XX/```:bangbang:
## 1. Load the env and Download the DB
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

#### 2-2. Change the format from fastq to fasta
```
source /programs/miniconda3/bin/activate run_dbcan4
export PATH=/programs/seqtk:$PATH

# need to do the below for all fastq, recommend to make a while loop as above linking process
seqtk seq -a input/501.fastq > input/501.fasta 
```

## 3. Run dbcan
:heavy_exclamation_mark: Methods the dbCAN uses are describe [here](https://dbcan.readthedocs.io/en/latest/user_guide/run_from_raw_reads.html)

#### 3-1. Use all methods 
```
run_dbcan input/501.fasta prok  --hmm_cpu 32 --out_dir /workdir/XXX/run_dbcan/output --use_signalP=TRUE -sp /programs/signalp-4.1/signalp --dia_cpu 32  --dbcan_thread 32 
```

#### 3-2. Use only 'HMMER vs dbCAN HMMdb' method 
```
run_dbcan input/501.fasta prok  --hmm_cpu 32 --out_dir /workdir/jp2626/run_dbcan/output --use_signalP=TRUE -sp /programs/signalp-4.1/signalp  --tools hmmer
```
