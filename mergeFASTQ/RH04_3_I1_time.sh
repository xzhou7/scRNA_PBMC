#!/bin/bash

echo "Start_combining_files"
date +"%T"

cat /Volumes/Xin_Lab/Single\ Cell\ Raw\ file/190726_COOPER_0285_AH2J2WBBXY_L1/RH04-3/Fastqs/*I1_001.fastq.gz /Volumes/Xin_Lab/Single\ Cell\ Raw\ file/190731_COOPER_0287_AH2GHFBBXY_L1/RH04-3/Fastqs/*I1_001.fastq.gz /Volumes/Xin_Lab/Single\ Cell\ Raw\ file/190731_COOPER_0287_AH2GHFBBXY_L2/RH04-3/Fastqs/*I1_001.fastq.gz /Volumes/Xin_Lab/Single\ Cell\ Raw\ file/190828_COOPER_0290_AH2HF5BBXY_L1/RH04-3/Fastqs/*I1_001.fastq.gz /Volumes/Xin_Lab/Single\ Cell\ Raw\ file/190828_COOPER_0290_AH2HF5BBXY_L2/RH04-3/Fastqs/*I1_001.fastq.gz /Volumes/Xin_Lab/Single\ Cell\ Raw\ file/190828_COOPER_0290_AH2HF5BBXY_L3/RH04-3/Fastqs/*I1_001.fastq.gz /Volumes/Xin_Lab/Single\ Cell\ Raw\ file/190828_COOPER_0290_AH2HF5BBXY_L4/RH04-3/Fastqs/*I1_001.fastq.gz /Volumes/Xin_Lab/Single\ Cell\ Raw\ file/190828_COOPER_0290_AH2HF5BBXY_L8/RH04-3/Fastqs/*I1_001.fastq.gz > /Volumes/Xin_Lab/SingleCell_GEO_PBMC/RH04-3/RH04-3_I1_001.fastq.gz

echo "file_moving_finished"
date +"%T"