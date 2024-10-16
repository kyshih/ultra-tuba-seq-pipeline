#!/usr/bin/env bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --mail-user=kyshih@stanford.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50g 
#SBATCH --time=24:00:00
#SBATCH --account=mwinslow
#SBATCH --partition=batch
#SBATCH --output=slurm_output/step5_%A.out

# activate conda env
module load adapterremoval/2.3.1
source /scg/apps/software/jupyter/python_3.9/etc/profile.d/conda.sh
conda activate UltraSeq

# Access environment variables
LP=${LP}
Input_experiment_ID=${Input_experiment_ID}
Project_directory=${Project_directory}
Python_script_address=${Python_script_address}
Step2_address=${Step2_address}
Step5_address=${Step5_address}

# Step5 
# Combined all the data
python3 $Python_script_address/UltraSeq_Step5.py --a "$Step2_address" \
--o "$Step5_address/"

echo "Step5. All the clustered gRNA and barcode data has been combined"
