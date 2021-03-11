#!/bin/bash

#sbatch Parameters
#SBATCH --time=4-00:00:00
#SBATCH --mail-user=pmaguire@stanford.edu
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_80
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=gssc
#SBATCH --account=org_gssc
#SBATCH --job-name RH03_Rebecca_Harper
#SBATCH --output Rebecca_Harper_Count_RH03_191017-13:59.out
#SBATCH --chdir /gssc/Users/pmaguire/Software/CellRanger_Sample_Processing/Logs
#SBATCH --mem=400G

#Loads Modules
module add legacy
module add bcl2fastq/2.20.0.422
module add cellranger/3.0.2

#Moves Into Base Folder
cd /gssc/Users/pmaguire/Rebecca_Harper/Results
mkdir -p RH03
cd RH03

#Count
echo "Sorting Fastq Files"

#190726_COOPER_0285_AH2J2WBBXY_L1 | RH03
cd /gssc/Users/pmaguire/Rebecca_Harper/Results/190726_COOPER_0285_AH2J2WBBXY_L1
mkdir -p RH03/Fastqs
mv RH03_*.fastq.gz RH03/Fastqs
#190731_COOPER_0287_AH2GHFBBXY_L1 | RH03
cd /gssc/Users/pmaguire/Rebecca_Harper/Results/190731_COOPER_0287_AH2GHFBBXY_L1
mkdir -p RH03/Fastqs
mv RH03_*.fastq.gz RH03/Fastqs
#190731_COOPER_0287_AH2GHFBBXY_L2 | RH03
cd /gssc/Users/pmaguire/Rebecca_Harper/Results/190731_COOPER_0287_AH2GHFBBXY_L2
mkdir -p RH03/Fastqs
mv RH03_*.fastq.gz RH03/Fastqs
#190828_COOPER_0290_AH2HF5BBXY_L1 | RH03
cd /gssc/Users/pmaguire/Rebecca_Harper/Results/190828_COOPER_0290_AH2HF5BBXY_L1
mkdir -p RH03/Fastqs
mv RH03_*.fastq.gz RH03/Fastqs
#190828_COOPER_0290_AH2HF5BBXY_L2 | RH03
cd /gssc/Users/pmaguire/Rebecca_Harper/Results/190828_COOPER_0290_AH2HF5BBXY_L2
mkdir -p RH03/Fastqs
mv RH03_*.fastq.gz RH03/Fastqs
#190828_COOPER_0290_AH2HF5BBXY_L3 | RH03
cd /gssc/Users/pmaguire/Rebecca_Harper/Results/190828_COOPER_0290_AH2HF5BBXY_L3
mkdir -p RH03/Fastqs
mv RH03_*.fastq.gz RH03/Fastqs
#190828_COOPER_0290_AH2HF5BBXY_L4 | RH03
cd /gssc/Users/pmaguire/Rebecca_Harper/Results/190828_COOPER_0290_AH2HF5BBXY_L4
mkdir -p RH03/Fastqs
mv RH03_*.fastq.gz RH03/Fastqs
#190828_COOPER_0290_AH2HF5BBXY_L8 | RH03
cd /gssc/Users/pmaguire/Rebecca_Harper/Results/190828_COOPER_0290_AH2HF5BBXY_L8
mkdir -p RH03/Fastqs
mv RH03_*.fastq.gz RH03/Fastqs

echo "Mapping Sample"

cd /gssc/Users/pmaguire/Rebecca_Harper/Results/RH03
cellranger count --id=RH03 --transcriptome=/gssc/Users/pmaguire/10X/Data/References/refdata-cellranger-GRCh38-1.2.0 --fastqs=/gssc/Users/pmaguire/Rebecca_Harper/Results/190726_COOPER_0285_AH2J2WBBXY_L1/RH03/Fastqs,/gssc/Users/pmaguire/Rebecca_Harper/Results/190731_COOPER_0287_AH2GHFBBXY_L1/RH03/Fastqs,/gssc/Users/pmaguire/Rebecca_Harper/Results/190731_COOPER_0287_AH2GHFBBXY_L2/RH03/Fastqs,/gssc/Users/pmaguire/Rebecca_Harper/Results/190828_COOPER_0290_AH2HF5BBXY_L1/RH03/Fastqs,/gssc/Users/pmaguire/Rebecca_Harper/Results/190828_COOPER_0290_AH2HF5BBXY_L2/RH03/Fastqs,/gssc/Users/pmaguire/Rebecca_Harper/Results/190828_COOPER_0290_AH2HF5BBXY_L3/RH03/Fastqs,/gssc/Users/pmaguire/Rebecca_Harper/Results/190828_COOPER_0290_AH2HF5BBXY_L4/RH03/Fastqs,/gssc/Users/pmaguire/Rebecca_Harper/Results/190828_COOPER_0290_AH2HF5BBXY_L8/RH03/Fastqs --sample=RH03 --expect-cells=5000  --localcores=24 --localmem=395 --chemistry=SC3Pv3

echo "Finished"
