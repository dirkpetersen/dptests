#!/bin/bash

while IFS= read -r line
do
    python3 eb_conda_configs/eb_conda_configs.py module --packages bioconda/"$line"
done < bioconda-packages.txt


