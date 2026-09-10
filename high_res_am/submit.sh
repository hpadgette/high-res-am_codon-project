#!/bin/bash

#SBATCH --qos=blanca-chbe-rdi
#SBATCH --time=12:00:00
#SBATCH --ntasks=4
#SBATCH --job-name=AM_tune
#SBATCH --constraint="cascadelake"
#SBATCH --output=final_tune.%j.out

source ~/.bashrc
conda activate cenv3_karplus
export LC_LANG=en_US.utf-8
export LANG=en_US.utf-8
echo "we're running AM!"
python -u -m run --config=config/config_example.json
