#!/bin/bash

pstacks -t sam -f ../batch1_bowtie2/PS12_018.sam -o ../batch1_bowtie2 -i 900 -m 10 -p 6 --model_type bounded 2>> ../batch1_bowtie2/pstacks_out_b1b2_wgenome

pstacks -t sam -f ../batch1_bowtie2/UP03_001.sam -o ../batch1_bowtie2 -i 900 -m 10 -p 6 --model_type bounded 2>> ../batch1_bowtie2/pstacks_out_b1b2_wgenome

pstacks -t sam -f ../batch1_bowtie2/PWS12_133.sam -o ../batch1_bowtie2 -i 900 -m 10 -p 6 --model_type bounded 2>> ../batch1_bowtie2/pstacks_out_b1b2_wgenome


sstacks -b 1 -c ../batch1_bowtie2/batch_1 -s ../batch1_bowtie2/PS12_018 -o ../batch1_bowtie2 -p 6 2>> sstacks_out_b1b2_wgenome
sstacks -b 1 -c ../batch1_bowtie2/batch_1 -s ../batch1_bowtie2/UP03_001 -o ../batch1_bowtie2 -p 6 2>> sstacks_out_b1b2_wgenome
sstacks -b 1 -c ../batch1_bowtie2/batch_1 -s ../batch1_bowtie2/PWS12_133 -o ../batch1_bowtie2 -p 6 2>> sstacks_out_b1b2_wgenome


