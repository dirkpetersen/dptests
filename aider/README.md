# Use Aider.chat with AWS Bedrock 

[Aider.chat](https://aider.chat/) is an AI based pair programmer that can write (almost) all your code for you. It works very well with [AWS Bedrock](https://aws.amazon.com/bedrock) and can run in the terminal Window of VS Code 

## Prepare AWS profile for Bedrock 

If you have not setup AWS yet, you need to get AWS credentials from your AWS Admin, ideally ask for a dedicated IAM user that only has AWSBedrockFullAccess permissions. Once you get your `access key` and `secret key` credentials for your IAM user create a new AWS profile just for Bedrock. 

Run the `aws` command in your Terminal on Linux/Mac or Powershell in Windows. If the `aws` command is not there, install it as [documented here](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) 

```
aws configure --profile mybedrock

AWS Access Key ID [None]: YOUR_ACCESS_KEY
AWS Secret Access Key [None]: YOUR_SECRET_KEY
Default region name [None]: us-west-2
Default output format [None]: 
```

If you work in an enterprise environment that requires [AWS SSO CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-sso.html), your AWS Admin may refuse to create a dedicated IAM user for Bedrock. In that case ask for Bedrock permissions AWSBedrockFullAccess to be added to your default setup (account/profile). 

Note: In addition to Bedrock permissions you or your AWS Admin also need to request access to each of the dozens of LLM models in your AWS Account. You only need to do this once for your team.
[Here is a 2 min Video on requesting model access](https://www.youtube.com/watch?v=WWHo7Awy0sQ)

## Install Aider and AWS boto3

Then we install Aider as well as the AWS API package boto3: 

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
# model: bedrock/us.anthropic.claude-sonnet-4-20250514-v1:0
# cache-prompts: true
```

Claude-3.5-Haiku is currently one of the fastest Anthropic models and is a good default, if you are an enterprise user with AWS SSO CLI and use your default profile, you can also put one of the most capable models (currently claude-sonnet-4) in the YML file and skip the next step

Now create a new custom .env config file "$HOME/.aider.aws.env" and add this content for best results in code writing:

```
AWS_PROFILE=mybedrock
AIDER_MODEL=bedrock/us.anthropic.claude-sonnet-4-20250514-v1:0
```

You can also replace the string `claude-sonnet-4` with `claude-opus-4` in the model name above to use an even more accurate model. However, Opus is more expensive and slower than Sonnet. In some other cases, you can achieve better results by splitting tasks between two different large language models: a reasoning model, such as Deepseek, responsible for the overall architecture, and a coding model that ensures correct code is written. In this case, we chose the previous Claude Sonnet 3.7 model because it is faster and slightly cheaper. Create a new custom .env config file "$HOME/.aider.aws-arch.env" and add this content 

```
AWS_PROFILE=mybedrock
AIDER_MODEL=bedrock/us.deepseek.r1-v1:0
AIDER_ARCHITECT=true
AIDER_EDITOR_MODEL=bedrock/us.anthropic.claude-3-7-sonnet-20250219-v1:0
```

Comment out or remove AWS_PROFILE if you are using AWS-SSO cli with the default AWS profile instead of a dedicated IAM user with an AWS profile mybedrock.  

## using Aider 

Now initialize a new git repository called `demo`:

```
mkdir demo; cd demo; git init
```

in which you launch aider (e.g. in the VS Code Terminal). If you use the `--env-file` on the CLI, you can easily switch between different models, for example by alternating between `"$HOME/.aider.aws-arch.env"` (reasoning for complex or custom architecture tasks) and `"$HOME/.aider.aws.env"` (claude only and faster). If you don't want to use the `--env-file` option, you can also edit a file called `.env` file at the root of your git repository, for example to better track which LLMs were used for this code. 

```
$ aider --env-file ~/.aider.aws.env
───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
You can skip this check with --no-gitignore
Add .aider* to .gitignore (recommended)? (Y)es/(N)o [Yes]: y
Added .aider* to .gitignore
Aider v0.83.2
Main model: bedrock/us.anthropic.claude-sonnet-4-20250514-v1:0 with diff edit format, infinite output
Weak model: bedrock/us.anthropic.claude-3-5-haiku-20241022-v1:0
Git repo: .git with 0 files
Repo-map: using 4096 tokens, auto refresh
Multiline mode: Enabled. Enter inserts newline, Alt-Enter submits text
───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
multi>

```

Now enter your prompt and send it to AWS Bedrock with ALT+Enter, for example :

`Make a very simple memory match game using python flask > v3 and modern web 2.0 tech`

Now open a new terminal, install Flask and run the Flask app

```
cd demo
python3 -m pip install -r requirements.txt
python3 app.py
```

Connect your web browser to `http://127.0.0.1:5000` then go back to the Aider terminal and ask: 

`Please change the colors of the tiles`

and confirm with `<ALT+Enter>`  

When you refresh your browser you will see that the color has changed. 


## more useful Aider options

you can also enter the model directly on the cli, for example if you are using the default AWS profile

```
aider --model bedrock/us.anthropic.claude-opus-4-20250514-v1:0
```

or if you prefer to interact with [Aider through a Browser](https://aider.chat/docs/usage/browser.html) instead of a terminal you can use the --browser option:

```
aider --env-file "$HOME/.aider.aws-r.env" --browser
```

For advanced use it is good practice to write the architecture including all requirements in a markdown file called CONVENTIONS.md and then launch that file with aider in read only mode to prevent that aider modifies this file. 


```
aider --env-file "$HOME/.aider.aws-r.env" --read CONVENTIONS.md
```

In the prompt state something like: `Please design the application as laid out in CONVENTIONS.md`. The benefits of this approach are documenting the architecture as fed into the Bedrock LLM and the ability to always remind the Bedrock LLM of the guidelines it must follow when coding.

IMPORTANT: If you enter anything in aider it will always modify code and add another git commit, however you can also enter `/undo <ALT+Enter>` to revert back to the previous git commit or you can use `/ask something <ALT+Enter>` to ask a question where no code will be changed. 

## addional Aider commands 

If you start your Aider terminal entry with a command such as /ask /undo /or run it will behave differently: 

- `/ask`: No code will be changed, Aider will just annwer your questions and give recommendations  
- `/undo`: Aider will undo the last git commmit (every single Aider code change will result in a git commit using a resonable git commit message)
- `/run`: will execute a command and parse the output into the Aider terminal  