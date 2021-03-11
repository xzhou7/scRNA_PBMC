#!/bin/bash

#sbatch Parameters
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=pmaguire@stanford.edu
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_80
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --cpus-per-task=12
#SBATCH --partition=gssc
#SBATCH --account=org_gssc
#SBATCH --job-name RH
#SBATCH --output Rebecca_Harper_Aggregation_RH_191022-18:13.out
#SBATCH --chdir /gssc/Users/pmaguire/Software/CellRanger_Sample_Processing/Logs
#SBATCH --mem=300G

#Loads Modules
module add legacy
module add bcl2fastq/2.20.0.422
module add cellranger/3.0.2

#Moves Into Base Folder
cd /gssc/Users/pmaguire/Rebecca_Harper/Results
mkdir -p RH
cd RH

#Aggregate
echo "Aggregating Samples"
cellranger aggr --id=RH --csv=/gssc/Users/pmaguire/Rebecca_Harper/Data/RH_aggregation_info.csv --normalize=mapped --localcores=12 --localmem=290

#Uploading Data
#Tarring Data
cd /gssc/Users/pmaguire/Rebecca_Harper/Results/RH
tar -czvf RH.tar.gz /gssc/Users/pmaguire/Rebecca_Harper/Results/RH/RH

#Uploading To Google Drive
module load gdrive/2.1
gdrive upload /gssc/Users/pmaguire/Rebecca_Harper/Results/RH/RH.tar.gz

echo "Finished"
