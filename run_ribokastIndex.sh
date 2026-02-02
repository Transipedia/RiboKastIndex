#!/bin/bash
#SBATCH --job-name="RiboKast_Index"
#SBATCH --partition=ssfa -t 100:00:00 --mem 100G
#SBATCH --cpus-per-task=6
source /home/safa.maddouri/miniconda3/bin/activate RiboKastIndex_env

SNAKEMAKE="snakemake"
SMK_DIR="/store/EQUIPES/SSFA/MEMBERS/safa.maddouri/RiboKastIndex_test"
SMK=${SMK_DIR}/snakefile
CONFIG=${SMK_DIR}/config.yaml

$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 3  --use-conda --unlock  
$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 3 --use-conda --rerun-incomplete --printshellcmds --reason

$SNAKEMAKE -s $SMK --configfile $CONFIG --dag | dot -Tpdf > /store/EQUIPES/SSFA/MEMBERS/safa.maddouri/RiboKastIndex_test/dag_pipeline.pdf

$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 3  --use-conda --reason 
