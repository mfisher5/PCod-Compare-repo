#!/bin/bash


bowtie2 -q -x ../refgenome/batch_7_ref_genome_b2 -U /media/mfisher5/New\ Volume/Mary/Stacks/PCod-US-repo/samplesT92/PS12_018.fq -S ../batch1_bowtie2/PS12_018.sam

bowtie2 -q -x ../refgenome/batch_7_ref_genome_b2 -U /media/mfisher5/New\ Volume/Mary/Stacks/PCod-US-repo/samplesT92/UP03_001.fq -S ../batch1_bowtie2/UP03_001.sam

bowtie2 -q -x ../refgenome/batch_7_ref_genome_b2 -U /media/mfisher5/New\ Volume/Mary/Stacks/PCod-US-repo/samplesT92/PWS12_133.fq -S ../batch1_bowtie2/PWS12_133.sam
