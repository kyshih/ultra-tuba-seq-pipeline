#!/usr/bin/env bash
#SBATCH --mail-user=kyshih@stanford.edu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=300G
#SBATCH --time=2:00:00
#SBATCH --account=mwinslow
#SBATCH --partition=batch
#SBATCH --output=slurm_output/step1_to_4_%A_%a.out

# Things to check:
# The fastq file address file, remember the text files must end with a newline
# sgRNA reference file
# Do you want to check the distance of sgRNA reference
# UltraSeq_Step3.py --a4 (4 is default)
# remember to change the array flag to match number of samples in the NGS_address
# make sure the position of the pair-end read files are correct. r1,r2,sampleID
# raw ultraseq data reside in NGS_Raw_data
# NGS address format is NGS_Raw_data/Input_experiment_ID/01.RawData/
# NGS address file is in Analysis_date

# modules
module load adapterremoval/2.3.1
#source /home/kyshih/.conda/pkgs/conda-4.13.0-py37h06a4308_0/etc/profile.d/conda.sh
source /scg/apps/software/jupyter/python_3.9/etc/profile.d/conda.sh
conda activate UltraSeq

# Access environment variables
LP=${LP}
Input_experiment_ID=${Input_experiment_ID}
Project_directory=${Project_directory}
Python_script_address=${Python_script_address}
Input_data_info_address=${Input_data_info_address}
Step1_address=${Step1_address}
Step2_address=${Step2_address}
Step5_address=${Step5_address}

echo $Python_script_address

# create directories for each step
mkdir -p $Step1_address
mkdir -p $Step2_address
mkdir -p $Step5_address

# Process each sample in parallel using job arrays
echo "SLURM_ARRAY_TASK_ID:" $SLURM_ARRAY_TASK_ID
sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p $Input_data_info_address)
# sample=$(sed -n '\$SLURM_ARRAY_TASK_ID\p' $Input_data_info_address)

# Step0: get read1, read2, and sampleID from NGS_address. each line is a sample aka mouse
r1=$(echo "$sample" | cut -d',' -f1)
r2=$(echo "$sample" | cut -d',' -f2)
sampleID=$(echo "$sample" | cut -d',' -f3)
echo "r1 ${r1}"
echo "r2 ${r2}"
echo "Step0. Retrieving read1 and read2 address completed for sample $sampleID"

# Step1: remove adaptors from both the 3' end of the reads
# each sample gets its own dir in the Merging dir
# basename specifies the perfix for the output files
temp_folder1=$Step1_address/$sampleID # Merging/sampleID
mkdir -p $temp_folder1
AdapterRemoval --file1 $r1  --file2 $r2 \
--adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
--adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT \
--basename "$temp_folder1/Merged"  --collapse --gzip
echo "Step1. For sample ${sampleID}, sequences merging is finished."

# Step2: generate bartender inputs. each sample gets its own dir in the Bartender dir
temp_folder2=$Step2_address/$sampleID # Bartender/sampleID
mkdir -p $temp_folder2
python3 $Python_script_address/UltraSeq_Step2.py --a "$temp_folder1/Merged.collapsed.gz" \
--o "$temp_folder2"
echo "Step2. gRNA bartender inputs were generated for ${sampleID}"

# Step3: cluster gRNAs using Bartender and generate Bartender_input_address file for 
#        each gRNA for clustering clonal barcodes in Step4
   # -f input file
   # -o output file prefix
   # -z merging threshold: -1 disable merging threshold to allow merging based on cluster distance only
   # -d max cluster length that may be merged. defult is 2. If the distance between two cluster sequences
   # is within this threshold, they will be merged 
   # -l seed length 5 as default
   # Bartender outputs 3 files: gRNA_cluster.csv, gRNA_barcode.csv, gRNA_quality.csv
   # see notion or webpage for details
bartender_single_com -z -1 -d 2 -l 5 -f "$temp_folder2/gRNA.bartender" \
-o "$temp_folder2/gRNA"
temp_folder3=$temp_folder2/Clonal_barcode
mkdir -p $temp_folder3

python3 $Python_script_address/UltraSeq_Step3.py --a1 "$temp_folder2/gRNA_barcode.csv" \
--a2 "$temp_folder2/gRNA_cluster.csv" --a3 "$Project_directory/gRNA_information.csv" \
--a5 "$temp_folder2/gRNA.bartender" --a6 "$temp_folder2/clonalbarcode.bartender" \
--o "$temp_folder2/"

echo "Step3. gRNA Clutering and generated Bartender_input_address for each gRNA has been completed"

# Step4: cluster clonal barcodes using Bartender
# Bartender_input_adress generated from Step3 has the address to all barcodes associated with a gRNA
# process all gRNAs in parallel
while read -r line2;
do 
   new_name=${line2/.bartender/}
   # echo $line2
   # echo $new_name
   bartender_single_com -z -1 -d 1 -l 5 -f "$line2" -o "$new_name"
done < $Step2_address/$sampleID/Bartender_input_address

echo "Step4. Barcode clustering has been completed for ${sampleID}"

# Step5
temp_folder4=$Step5_address/$sampleID
mkdir -p $temp_folder4


