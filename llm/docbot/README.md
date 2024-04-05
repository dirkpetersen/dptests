
# Run Llama2 7BN model over some wiki docs to generate a chatbot

use doc: https://swharden.com/blog/2023-07-30-ai-document-qa/

wget https://huggingface.co/TheBloke/Llama-2-7B-Chat-GGML/resolve/main/llama-2-7b-chat.ggmlv3.q8_0.bin
mkdir -p ./wiki-acc
wget -r -l10 -np -e robots=off -P ./wiki-acc https://wiki.ohsu.edu/display/ACC


# Run a service on slurm 

Check out the slurm service controller [ssvc](../../ssvc/ssvc.md) 

