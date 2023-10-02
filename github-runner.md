# installing a github runner 

## Repos setup

on the machine that hosts the runner execute: 

```
mkdir -p ~/gh && cd ~/gh
export RUNNER_ALLOW_RUNASROOT=1
```

* In your githus repos go to settings -> Actions -> Runners and setup a new self hosted runner and then paste all the commands into your shell 

* this will be something similar to 

```
mkdir actions-runner && cd actions-runner
curl -o actions-runner-linux-x64-2.309.0.tar.gz -L https://github.com/actions/runner/releases/download/v2.309.0/actions-runner-linux-x64-2.309.0.tar.gz
tar xzf ./actions-runner-linux-x64-2.309.0.tar.gz
./config.sh --url https://github.com/dirkpetersen/dptests --token XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
```

* Do NOT execute run.sh instead install and run it as a system service 

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