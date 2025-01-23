from flask import Flask, render_template, request
import boto3
import PyPDF2
import io
import os
import json
from werkzeug.utils import secure_filename
from langchain_community.embeddings import BedrockEmbeddings
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.vectorstores import FAISS

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max file size
app.config['SECRET_KEY'] = os.urandom(24)

# Configure AWS
bedrock = boto3.client(
    service_name='bedrock-runtime',
    region_name='us-east-1'
)

# Configure embeddings and text splitter
embeddings = BedrockEmbeddings(
    client=bedrock,
    model_id="amazon.titan-embed-text-v1"
)

text_splitter = RecursiveCharacterTextSplitter(
    chunk_size=1000,
    chunk_overlap=100
)

def extract_text_from_pdf(pdf_file):
    pdf_reader = PyPDF2.PdfReader(pdf_file)
    text = ""
    for page in pdf_reader.pages:
        text += page.extract_text()
    return text

def evaluate_requirements(policy_text, submission_text):
    # Split documents into chunks
    policy_chunks = text_splitter.split_text(policy_text)
    submission_chunks = text_splitter.split_text(submission_text)
    
    # Create vector store from policy document
    policy_store = FAISS.from_texts(policy_chunks, embeddings)
    
    # For each submission chunk, find relevant policy requirements
    relevant_pairs = []
    for chunk in submission_chunks:
        similar_docs = policy_store.similarity_search(chunk, k=2)
        relevant_pairs.extend([(chunk, doc.page_content) for doc in similar_docs])
    
    # Create analysis prompt with relevant chunks
    analysis_prompt = """Human: I will provide pairs of text chunks from two documents. 
For each pair, the first is from a submission document and the second is from a policy document. 
Analyze if the submission meets the requirements in the policy.\n\n"""
    
    for sub_chunk, pol_chunk in relevant_pairs:
        analysis_prompt += f"Submission chunk:\n{sub_chunk}\n\nMatching policy chunk:\n{pol_chunk}\n\n"
    
    analysis_prompt += "Based on all these comparisons, respond with exactly one word (GREEN, YELLOW, or RED) "
    analysis_prompt += "followed by a brief explanation if YELLOW. GREEN means the submission fully meets requirements, "
    analysis_prompt += "RED means it clearly doesn't, and YELLOW means there are uncertainties that need human review."
    
    response = bedrock.invoke_model(
        modelId="anthropic.claude-v3-sonnet",
        body=json.dumps({
            "anthropic_version": "bedrock-2023-05-31",
            "max_tokens": 1000,
            "temperature": 0,
            "messages": [
                {
                    "role": "user",
                    "content": analysis_prompt
                }
            ]
        })
    )
    
    result = json.loads(response['body'].read())['completion']
    
    if "GREEN" in result.upper():
        return "GREEN", ""
    elif "RED" in result.upper():
        return "RED", ""
    else:
        # Extract explanation after "YELLOW"
        parts = result.split("YELLOW", 1)
        explanation = parts[1].strip() if len(parts) > 1 else "No specific explanation provided"
        return "YELLOW", explanation

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        if 'policy' not in request.files or 'submission' not in request.files:
            return render_template('index.html', error="Both files are required")
        
        policy_file = request.files['policy']
        submission_file = request.files['submission']
        
        if policy_file.filename == '' or submission_file.filename == '':
            return render_template('index.html', error="Both files are required")
        
        try:
            policy_text = extract_text_from_pdf(policy_file)
            submission_text = extract_text_from_pdf(submission_file)
            
            result, explanation = evaluate_requirements(policy_text, submission_text)
            return render_template('index.html', result=result, explanation=explanation)
            
        except Exception as e:
            return render_template('index.html', error=str(e))
    
    return render_template('index.html')

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
