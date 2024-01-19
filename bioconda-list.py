#! /usr/bin/env python3
import subprocess

class BigBadClass:
    @staticmethod
    def execute_and_parse(cmd):
        # Run the command and capture the output
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()

        # Check for errors
        if process.returncode != 0:
            raise Exception(f"Error executing command: {stderr.decode('utf-8')}")

        # Parse the output into tuples
        lines = stdout.decode('utf-8').strip().split('\n')
        tuples = [tuple(line.split()) for line in lines if line]

        return tuples

# Example usage
cmd = "pixi --color never search -c bioconda -q -l 25000 '*'"
big_bad = BigBadClass()
result = big_bad.execute_and_parse(cmd)
for line in result:
    print(line[0])
#print(result)
