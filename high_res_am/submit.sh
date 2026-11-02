#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --qos=blanca-chbe-rdi
#SBATCH --account=blanca-chbe-rdi
#SBATCH --constraint="cascadelake"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --job-name=gl-AM-test3
#SBATCH --output=CLK19gl_eOD-GT8_AM_test3.%j.out
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=hapa3524@colorado.edu

echo this is job $SLURM_JOBID
echo $SLURM_NTASKS

source ~/kayla-conda.sh
conda activate cenv3_karplus

export LC_LANG=en_US.utf-8
export LANG=en_US.utf-8

echo "we're running AM!"
python -u -m run --config=codon-project/config_CLK19gl.json