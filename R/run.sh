#! /bin/bash

DEBUG=0

ml fhR/4
echo using R version:
which R

mkdir -p $2
for f in $(ls -1 $1); do
  base=$(basename $f)
  if [[ ${DEBUG} -eq 0 ]]; then
    sbatch -p restart-new --requeue -c 1 think.R $1/${f} $2/out-${base}
  else
    echo "input: $1/${f}, output: $2/out-${base}"
    ./think.R $1/${f} $2/out-${base}
  fi
done

