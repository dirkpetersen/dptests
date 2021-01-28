#! /usr/bin/env Rscript

restartcount <- as.integer(Sys.getenv('SLURM_RESTART_COUNT', 0))
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]

if (restartcount > == 0) {
  if (file.exists(outfile)) {
    file.remove(outfile)
  }
}

df1 <- read.csv(file = infile)
write.csv(df1, outfile, row.names = FALSE)
