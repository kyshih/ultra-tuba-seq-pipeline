#!/usr/bin/env bash
#SBATCH --mail-user=kyshih@stanford.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50g
#SBATCH --time=24:00:00
#SBATCH --account=mwinslow
#SBATCH --partition=batch
#SBATCH --output=slurm_output/ultraseq_pipeline-%j.out

# Overall directory and project info
LP="/labs/mwinslow/Karen/Ultra_seq/"
Input_experiment_ID="test"
Project_directory=$LP$Input_experiment_ID
Python_script_address=$LP$Input_experiment_ID

# Input and output directories
Input_data_info_address="$Project_directory/NGS_address"
Step1_address="$Project_directory/Merging"
Step2_address="$Project_directory/Bartender"
Step5_address="$Project_directory/Processed_data"

number_of_samples=10

# Export environment variables to be accessible in Step 1-4 and Step 5 scripts
export LP Input_experiment_ID Project_directory Python_script_address Input_data_info_address Step1_address Step2_address Step5_address

# Activate conda environment (if needed)
module load adapterremoval/2.3.1
source /scg/apps/software/jupyter/python_3.9/etc/profile.d/conda.sh
conda activate UltraSeq

# Submit the job array for Steps 1-4 (each sample), passing the necessary arguments
jid1=$(sbatch --array=1-$number_of_samples --parsable UltraSeq_Step1_to_4.sh)
echo "Submitted Step 1-4 job array with Job ID $jid1"

# Submit the Step 5 script with dependency on the completion of all Step 1-4 jobs
jid2=$(sbatch --dependency=afterok:$jid1 --parsable UltraSeq_Step5.sh)
#jid2=$(sbatch --parsable UltraSeq_Step5.sh)
echo "Submitted Step 5 job ID: $jid2, dependent on Step 1-4 completion"
