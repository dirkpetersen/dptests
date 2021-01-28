## using the restart queue 

This is an example for cleanly using the restart queue. 

Start the script with this command:

```
./run.sh ./csv ./out
```

The script will :

* loop through folder csv
* lauch script think.R with each data file from the csv folder and a new output csv file as 2nd argument
* if `SLURM_RESTART_COUNT` is > 0  ensure that the output CSV file is deleted to allow a re-queued job to start from scratch* the output csv files should contain only the first few lines of the original csv files




