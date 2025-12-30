#!/bin/bash
# MODIFIED FROM https://github.com/RyanLeeHayes/ALF/blob/master/examples/engines/blade/subsetAll.sh

start_time=$(date +%s)
echo "Script started at $(date)"

step1_start=$(date +%s)

export step=112
echo '>>> running  5ns production simulation with 5 trials(a,b,c,d,e)'
for p in a b c d e # N=5
do
  echo ' >> running letter' $p
  export p
  export nitt=5 # was 5
  #PID=`sbatch --parsable $SLURMOPTSMD --array=1-1%1 $DEPEND ./runprod.sh`
  #NEWDEPEND="${NEWDEPEND}${COMMA}afterok:$PID"
  #COMMA=","
  for ii in {1..1} # So S=5ns
  do
    echo '  > running ii=' $ii
    export SLURM_ARRAY_TASK_ID=$ii
    iend=$(( $SLURM_ARRAY_TASK_ID * $nitt ))
    ibeg=$(( $iend - $nitt ))
    echo '  > ibeg:' $ibeg 'iend:' $iend
    python -c "import alf; alf.runprod($step,'$p',$ibeg,$iend,nsteps=500000,engine='blade')"
  done
done
step1_end=$(date +%s)
step1_elapsed=$((step1_end - step1_start))
step1_hours=$((step1_elapsed / 3600))
step1_minutes=$((step1_elapsed % 3600 / 60))
echo "[TIMING] this step took $step1_hours hours $step1_minutes minutes"


step1_start=$(date +%s)

echo '>>> postprocessing...'
export i=112
export eqS=1
export S=5
export N=5
export skipE=1
python -c "import alf; alf.postprocess($i,$eqS,$S,$N,$skipE,True,engine='blade',G_imp='../prep/G_imp',lc=0.0)"

step1_end=$(date +%s)
step1_elapsed=$((step1_end - step1_start))
step1_hours=$((step1_elapsed / 3600))
step1_minutes=$((step1_elapsed % 3600 / 60))
echo "[TIMING] this step took $step1_hours hours $step1_minutes minutes"


finish_time=$(date +%s)
echo "Script finished at $(date)"


