#!/bin/bash

chromosomes=(1 2 3 4 5 Mt Pt)
contexts=(CpG CHG CHH)

input_dir="temperatureComp/testSampleCSVs"

for context in "${contexts[@]}"; do
    context_dir="${input_dir}/${context}"
    
    if [ ! -d "$context_dir" ]; then
        mkdir -p "$context_dir"
    fi

    input_file="${input_dir}/testSamples${context^^}.csv"

    for chrom in "${chromosomes[@]}"; do
        output_file="${context_dir}/output_${context}_${chrom}.csv"

        if [ ! -f "$output_file" ]; then
            ./shannon div -m "$input_file" -o "$output_file" -s "$chrom" -c 5 6 -n m um
        else
            echo "File $output_file exists!"
        fi
    done
done

echo "All files are created!"
