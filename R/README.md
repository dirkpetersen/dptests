## using the restart queue 

This is an example for cleanly using the restart queue. 

Start the script with this command:

```
./run.sh ./csv ./out
```

The script will :

- loop through folder csv
- lauch script think.R with each data file from the csv folder and a new output csv file as 2nd argument
- if an interupted job is restarted (`SLURM_RESTART_COUNT` > 0) ensure that the incomplete output CSV file is deleted to allow a re-queued job to start from scratch
- in the example the output csv files should contain only the first few lines of the original csv files

