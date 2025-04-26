# Use Aider.chat with AWS Bedrock 

## Prepare AWS profile for Bedrock 

get AWS credentials from your AWS Admin, ideally for an IAM user that only has AWSBedrockFullAccess permissions. Then create a new AWS profile just for Bedrock. 

Run the `aws` command in your Terminal on Linux/Mac or Powershell in Windows. If the `aws` command is not there, install it as [documented here](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) 

```
aws configure --profile mybedrock

AWS Access Key ID [None]: YOUR_ACCESS_KEY
AWS Secret Access Key [None]: YOUR_SECRET_KEY
Default region name [None]: us-west-2
Default output format [None]: 
```

If you work in an enterprise environment that requires [AWS SSO CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-sso.html), your AWS Admin may refuse to create an IAM user for you. In that case ask for Bedrock to be added to your default profile. 

## Install Aider and AWS boto3

Linux / Mac: 

```
curl -LsSf https://aider.chat/install.sh | sh
uv tool run --from aider-chat pip install boto3
```

Windows:

```
powershell -ExecutionPolicy ByPass -c "irm https://aider.chat/install.ps1 | iex"
uv tool run --from aider-chat pip install boto3
```

## configure Aider 

Create an empty global config file $HOME/.aider.conf.yml or download the template: 


Linux/Mac:

```
curl -o $HOME/.aider.conf.yml https://raw.githubusercontent.com/Aider-AI/aider/refs/heads/main/aider/website/assets/sample.aider.conf.yml
```

Windows:

```
iwr -Uri "https://raw.githubusercontent.com/Aider-AI/aider/refs/heads/main/aider/website/assets/sample.aider.conf.yml" -OutFile "$HOME/.aider.conf.yml"
```

Add this to the top of file "$HOME/.aider.conf.yml":

```
multiline: true
vim: true
dark-mode: true
show-model-warnings: false
subtree-only: true
model: bedrock/us.anthropic.claude-3-5-haiku-20241022-v1:0
# cache-prompts: true
```

Claude-3.5-Haiku is currently the fastest Anthropic model and is a good default 

create an .env config file "$HOME/.aider.aws-r.env" and add this content for best results in code writing:

```
AWS_PROFILE=mybedrock
AIDER_ARCHITECT=true
AIDER_MODEL=bedrock/us.deepseek.r1-v1:0
AIDER_EDITOR_MODEL=bedrock/us.anthropic.claude-3-7-sonnet-20250219-v1:0
```

for complex tasks that require reasoning ..... but for simpler tasks create another .env config file "$HOME/.aider.aws-c.env" and add this content: 

```
AWS_PROFILE=mybedrock
AIDER_ARCHITECT=false
AIDER_MODEL=bedrock/us.anthropic.claude-3-7-sonnet-20250219-v1:0
```

Comment out or remove AWS_PROFILE if you are using AWS-SSO cli with the default AWS profile   

Now initialize a new git repository called `demo`:

```
mkdir demo; cd demo; git init
```

in which you launch aider (e.g. in the VS Code Terminal). If you use the `--env-file` on the CLI, you can easily switch between different models, for example by alternating between `"$HOME/.aider.aws-r.env"` (reasoning) and `"$HOME/.aider.aws-c.env"` (claude). If you don't want to use the `--env-file` option, you can also edit a file called `.env` file at the root of your git repository, for example to better track which LLMs were used for this code. 

```
aider --env-file "$HOME/.aider.aws-r.env"

───────────────────────────────────────────────────────────────────────────────────────────────
Aider v0.82.2
Model: bedrock/us.deepseek.r1-v1:0 with architect edit format
Editor model: bedrock/us.anthropic.claude-3-7-sonnet-20250219-v1:0 with editor-diff edit format
Git repo: .git with 1,748 files
Repo-map: using 4096 tokens, files refresh
Multiline mode: Enabled. Enter inserts newline, Alt-Enter submits text
Cost estimates may be inaccurate when using streaming and caching.
───────────────────────────────────────────────────────────────────────────────────────────────
architect multi>

```

Now enter your prompt and send it to AWS Bedrock with ALT+Enter, for example :

`Make a simple memory match game using python flask and modern web 2.0 tech`

Now open a new terminal, install Flask and run the Flask app

```
cd demo
python3 -m pip install flask
python3 app.py
```

Connect your web browser to `http://127.0.0.1:5000` then go back to the Aider terminal and ask: 

`Please change the colors of the tiles`

and confirm with `<ALT+Enter>`

you can also enter the model directly on the cli, for example if you are using the default AWS profile

```
aider --model bedrock/us.anthropic.claude-3-7-sonnet-20250219-v1:0
```

For advanced use it is good practice to write the architecture including all requirements in a markdown file called CONVENTIONS.md and then launch that file with aider in read only mode to prevent that aider modifies this file. 


```
aider --env-file "$HOME/.aider.aws-r.env" --read CONVENTIONS.md
```

In the prompt state something like: `Please design the application as laid out in CONVENTIONS.md`. The benefits of this approach are documenting the architecture as fed into the Bedrock LLM and the ability to always remind the Bedrock LLM of the guidelines it must follow when coding.

IMPORTANT: If you enter anything in aider it will always modify code and add another git commit, however you can also enter `/undo <ALT+Enter>` to revert back to the previous git commit or you can use `/ask something <ALT+Enter>` to ask a question where no code will be changed. 

