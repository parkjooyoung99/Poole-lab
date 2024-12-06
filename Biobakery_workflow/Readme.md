# :bangbang: All work should be done in ```/workdir/XX/```:bangbang:
## 1. Set up env
#### 1-1. Download [env_biobakery_pip.yml]([Biobakery_workflow/env_biobakery_pip.yml](https://github.com/parkjooyoung99/Poole-lab/blob/33dfbbb9f957c28d6942f76c3527b276a539b8b5/Biobakery_workflow/env_biobakery_pip.yml)) file
#### 1-2. Create conda env in server
```
source $HOME/miniconda3/bin/activate
conda create -f pathtoyml/env_biobakery_pip.yml
conda activate biobakery_pip
```
#### 1-3. Download databases for wgmx ####
```
mkdir /workdir/jp2626/biobakery_workflow_databases
biobakery_workflows_databases --install wmgx --location  /workdir/jp2626/biobakery_workflow_databases
export BIOBAKERY_WORKFLOWS_DATABASES=/workdir/jp2626/biobakery_workflow_databases
```
## 2. Run your data
:heavy_exclamation_mark: This is for wgmx. If your data isn't wgmx check [biobakery_workflow github](https://github.com/biobakery/biobakery_workflows)

#### 2-1. Rename your fastq.gz ####
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


