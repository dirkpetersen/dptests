# installing a github runner 

## Repos setup

on the machine that hosts the runner execute: 

```
mkdir -p ~/gh && cd ~/gh
export RUNNER_ALLOW_RUNASROOT=1
```

* In your githb repos go to settings -> Actions -> Runners and setup a new self hosted runner and then paste all the commands into your shell 

* Do NOT execute run.sh instead install and run it as a sysyemd service 

```
./svc.sh install
./svc.sh start
```

more info 

https://docs.github.com/en/actions/hosting-your-own-runners/managing-self-hosted-runners/configuring-the-self-hosted-runner-application-as-a-service



Using your self-hosted runner
# Use this YAML in your workflow file for each job
runs-on: self-hosted

## Github Actions 

* Go to Actions -> New workflow -> Skip this and set up a workflow yourself 
* paste the below into a new workflow main.yml 


```
name: Pull On Commit

on:
  push:
    branches:
      - main # or the name of your default branch

jobs:
  pull:
    runs-on: self-hosted
    steps:
      - name: Checkout Code
        uses: actions/checkout@v2
      
      - name: Pull Latest Code
        run: git pull origin main # or the name of your default branch

```