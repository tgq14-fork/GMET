#!/bin/bash

#! Define configuration flags
# See: https://www.acrc.bris.ac.uk/protected/bc4-docs/scheduler/index.html

#SBATCH --job-name=Clen
#SBATCH --time=0-1:0:0
#SBATCH --mem=20G

#! add the MATLAB module (as per BCp4)
module load matlab/R2017b
#! Run the job
matlab -nodisplay -r gmet_process_climo_ensemble