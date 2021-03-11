#!/bin/bash

#sbatch Parameters
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=pmaguire@stanford.edu
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_80
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --account=mpsnyder
#SBATCH --job-name 190828_COOPER_0290_AH2HF5BBXY_L2_Rebecca_Harper
#SBATCH --output Rebecca_Harper_Mkfastq_190828_COOPER_0290_AH2HF5BBXY_L2_191001-22:58.out
#SBATCH --chdir /home/pmaguire/Software/CellRanger_Sample_Processing/Logs
#SBATCH --mem=300G

#Loads Modules
module add legacy
module add bcl2fastq/2.20.0.422
module add cellranger/3.0.2

#Moves Into Base Folder
cd /seqctr/archive/gsfs0/seq_center/Users/pmaguire/Rebecca_Harper/Results
mkdir -p 190828_COOPER_0290_AH2HF5BBXY_L2
cd 190828_COOPER_0290_AH2HF5BBXY_L2

#Mkfastq
echo "Demultiplexing Samples"
cellranger mkfastq --run=/seqctr/archive/gsfs0/seq_center/Users/pmaguire/Rebecca_Harper/Data/190828_COOPER_0290_AH2HF5BBXY_L2 --csv=/seqctr/archive/gsfs0/seq_center/Users/pmaguire/Rebecca_Harper/Data/190828_COOPER_0290_AH2HF5BBXY_L2/190828_COOPER_0290_AH2HF5BBXY_L2_sample_sheet.csv --localcores=12 --localmem=300 --output-dir=/seqctr/archive/gsfs0/seq_center/Users/pmaguire/Rebecca_Harper/Results/190828_COOPER_0290_AH2HF5BBXY_L2

echo "Finished"
