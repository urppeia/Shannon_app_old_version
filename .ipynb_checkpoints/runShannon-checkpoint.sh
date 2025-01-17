#!/usr/bin/env python3
# -*- coding: utf-8 -*-

file_names = [
    "metadata_CHG.Col-0.wt.aerial-part.csv",
    "metadata_CHG.Col-0.wt.embryo.csv",
    "metadata_CHG.Col-0.wt.endosperm.csv",
    "metadata_CHG.Col-0.wt.immature-flower-buds.csv",
    "metadata_CHG.Col-0.wt.inflorescence.csv",
    "metadata_CHG.Col-0.wt.root.csv",
    "metadata_CHG.Col-0.wt.rosette.csv",
    "metadata_CHG.Col-0.wt.shoot.csv",
    "metadata_CHG.Col-0.wt.sperm-cell.csv",
    "metadata_CHG.Col-0.wt.vegetative-nucleus.csv",
    "metadata_CHG.Col-0.wt.whole-organism.csv",
    "metadata_CHH.Col-0.wt.aerial-part.csv",
    "metadata_CHH.Col-0.wt.embryo.csv",
    "metadata_CHH.Col-0.wt.endosperm.csv",
    "metadata_CHH.Col-0.wt.immature-flower-buds.csv",
    "metadata_CHH.Col-0.wt.inflorescence.csv",
    "metadata_CHH.Col-0.wt.root.csv",
    "metadata_CHH.Col-0.wt.rosette.csv",
    "metadata_CHH.Col-0.wt.shoot.csv",
    "metadata_CHH.Col-0.wt.sperm-cell.csv",
    "metadata_CHH.Col-0.wt.vegetative-nucleus.csv",
    "metadata_CHH.Col-0.wt.whole-organism.csv",
    "metadata_CpG.Col-0.wt.aerial-part.csv",
    "metadata_CpG.Col-0.wt.embryo.csv",
    "metadata_CpG.Col-0.wt.endosperm.csv",
    "metadata_CpG.Col-0.wt.immature-flower-buds.csv",
    "metadata_CpG.Col-0.wt.inflorescence.csv",
    "metadata_CpG.Col-0.wt.root.csv",
    "metadata_CpG.Col-0.wt.rosette.csv",
    "metadata_CpG.Col-0.wt.shoot.csv",
    "metadata_CpG.Col-0.wt.sperm-cell.csv",
    "metadata_CpG.Col-0.wt.vegetative-nucleus.csv",
    "metadata_CpG.Col-0.wt.whole-organism.csv",
    "metadata_handcurated.csv",
    "metadata_handcurated_marc.csv"
]


if [ ! -d "ShannonData" ]; then
    mkdir "ShannonData"
fi

for input_file in "${file_names[@]}"; do
    output_file="ShannonData/output_${input_file}" 
    sequence=1

    if [ ! -f "$output_file" ]; then
        ./shannon div -m "../meta-methylome/results/tables/metadata/${input_file}" -o "${output_file}" -s "$sequence" -c 4 5 6 -n cov m um 
    else
        echo "File $output_file exists!"
    fi
done

echo "All Files are created!"
