#!/bin/bash

files=(*.fasta)
first_100=("${files[@]::100}")

let i=1
total="${#first_100[@]}"

for file in "${first_100[@]}"
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

# took ~4 minutes for 100 alignments
