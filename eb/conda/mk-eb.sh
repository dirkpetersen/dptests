#!/bin/bash


/usr/bin/python3 -m pip install anaconda-client jinja2

while IFS= read -r line
do
    /usr/bin/python3 eb_conda_configs.py module --packages bioconda/"$line"
done < bioconda-packages-working-x86-64.txt

