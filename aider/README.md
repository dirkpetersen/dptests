# Use Aider.chat with AWS Bedrock 

## Prepare AWS profile for Bedrock 

get AWS credentials from your AWS Admin, 
ideally for an IAM user that only has AWSBedrockFullAccess permissions 
create a new AWS profile just for Bedrock

```
aws configure --profile mybedrock

AWS Access Key ID [None]: YOUR_ACCESS_KEY
AWS Secret Access Key [None]: YOUR_SECRET_KEY
Default region name [None]: us-west-2
Default output format [None]: 
```

## Install Aider

```
curl -LsSf https://aider.chat/install.sh | sh
```

## configure Aider 

Create an empty global config file ~/.aider.conf.yaml or download the template: 

```
curl -o ~/.aider.conf.yml https://raw.githubusercontent.com/Aider-AI/aider/refs/heads/main/aider/website/assets/sample.aider.conf.yml
```

Add this to the top of file ~/.aider.conf.yaml

```
multiline: true
dark-mode: true
cache-prompts: true
show-model-warnings: false
vim: true
```

create an .env config file ~/.aider.aws.env and add this content:

```
AWS_PROFILE=mybedrock
AIDER_ARCHITECT=true
AIDER_MODEL=bedrock/us.deepseek.r1-v1:0
AIDER_EDITOR_MODEL=bedrock/us.anthropic.claude-3-7-sonnet-20250219-v1:0
```

In the future you will end up using multiple .env files with different configurations.

Now enter a git repository in which you launch aider (e.g. in the VS Code Terminal). If you don't want to 
enter --env-file on the CLI you can also edit an '.env' file at the root of your git repository to track 
which LLMs were used for this code. 

```
aider --env-file ~/.aider.aws.env

───────────────────────────────────────────────────────────────────────────────────────────────
Aider v0.82.2
Model: bedrock/us.deepseek.r1-v1:0 with architect edit format
Editor model: bedrock/us.anthropic.claude-3-7-sonnet-20250219-v1:0 with editor-diff edit format
Git repo: .git with 1,748 files
Warning: For large repos, consider using --subtree-only and .aiderignore
See: https://aider.chat/docs/faq.html#can-i-use-aider-in-a-large-mono-repo
Repo-map: using 4096 tokens, files refresh
Multiline mode: Enabled. Enter inserts newline, Alt-Enter submits text
Cost estimates may be inaccurate when using streaming and caching.
───────────────────────────────────────────────────────────────────────────────────────────────
architect multi>

```

Now enter your prompt and send it to AWS Bedrock with ALT+Enter





