#!/bin/sh
#PBS -l select=1:ncpus=8:mem=500gb
#PBS -l walltime=150:00:00
#export TMPDIR=/data/work/I2BC/safa.maddouri/tmp
#cd /data/work/I2BC/safa.maddouri
#source /data/work/I2BC/safa.maddouri/miniconda3/bin/activate 

source /home/safa.maddouri/miniconda3/bin/activate Ribodocs_env

SNAKEMAKE="snakemake"
SMK_DIR="/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/RiboFramedKmers"
SMK=${SMK_DIR}/snakefile
CONFIG=${SMK_DIR}/config.yaml

$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 2  --use-conda --unlock  
$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 2  --use-conda --reason 
