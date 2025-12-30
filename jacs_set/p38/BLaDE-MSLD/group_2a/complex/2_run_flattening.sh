#!/bin/bash
# MODIFIED FROM https://github.com/RyanLeeHayes/ALF/blob/master/examples/engines/blade/subsetAll.sh

start_time=$(date +%s)
echo "Script started at $(date)"

PATCH_FILE="alf_postprocess.patch"
ALF_FILE=$(python -c "import alf; print(alf.postprocess.__code__.co_filename)")
if patch --dry-run --silent "$ALF_FILE" "$PATCH_FILE" > /dev/null 2>&1; then
    patch "$ALF_FILE" "$PATCH_FILE"
fi

python -c "import alf; alf.initialize(engine='blade',minimize=False)"


step1_start=$(date +%s)
python -c "import alf; alf.runflat(1,100,13000,39000,engine='blade',G_imp='../prep/G_imp')"
step1_end=$(date +%s)
step1_elapsed=$((step1_end - step1_start))
step1_hours=$((step1_elapsed / 3600))
step1_minutes=$((step1_elapsed % 3600 / 60))
echo "[TIMING] this step took $step1_hours hours $step1_minutes minutes"


step1_start=$(date +%s)
python -c "import alf; alf.runflat(101,111,125000,375000,engine='blade',G_imp='../prep/G_imp')"
step1_end=$(date +%s)
step1_elapsed=$((step1_end - step1_start))
step1_hours=$((step1_elapsed / 3600))
step1_minutes=$((step1_elapsed % 3600 / 60))
echo "[TIMING] this step took $step1_hours hours $step1_minutes minutes"


finish_time=$(date +%s)
echo "Script finished at $(date)"


