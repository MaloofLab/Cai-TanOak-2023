#!/bin/bash

files=(*.fasta)
first_1000=("${files[@]::1000}")

let i=1
total="${#first_1000[@]}"

for file in "${first_1000[@]}"
do
    echo "Aligning: $i/$total"

    out_nt="../aligned_seqs/${file}.out_nt"
    out_aa="../aligned_seqs/${file}.out_aa"
    java -jar ~/macse_v2.05.jar -prog alignSequences -seq $file -out_NT $out_nt -out_AA $out_aa

    #echo "$file"
    #echo "$out_nt"
    #echo "$out_aa"
    let i++
done

# 1000 jobs took ~48 minutes
