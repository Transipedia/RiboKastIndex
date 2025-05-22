set -e

# Involved programs
KAMRAT_IMG=$1 # need to change the singularity image path accordingly
TIME='/usr/bin/time -f "Command - %C\n\tuser time in seconds -  %U\n\tsystem (kernel) time in seconds - %S\n\telapsed real time (wall clock) in seconds - %e\n\tpercent of CPU this job got - %P"'


# Inputs
IN_MAT=$2

# Outputs
OUT_DIR=$3
IDX_DIR="$OUT_DIR""/index/"

kmerSize=$4
kmerOverlap=$((kmerSize-1))
mkdir -p $IDX_DIR
mkdir -p $OUT_DIR

bash -c "$TIME apptainer exec -B '/store:/store' -B '/data:/data' $KAMRAT_IMG kamrat index -intab $IN_MAT -outdir $IDX_DIR -klen $kmerSize "

bash -c "$TIME apptainer exec -B '/store:/store' -B '/data:/data' $KAMRAT_IMG kamrat merge -idxdir $IDX_DIR -overlap $kmerOverlap-15 -interv pearson:0.20 -withcounts mean -outpath $OUT_DIR/merged-res.tsv"
