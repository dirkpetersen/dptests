from flask import Flask, render_template, request
import boto3
import PyPDF2
import io
import os
from werkzeug.utils import secure_filename

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max file size
app.config['SECRET_KEY'] = os.urandom(24)

# Configure AWS
bedrock = boto3.client(
    service_name='bedrock-runtime',
    region_name='us-east-1'
)

def extract_text_from_pdf(pdf_file):
    pdf_reader = PyPDF2.PdfReader(pdf_file)
    text = ""
    for page in pdf_reader.pages:
        text += page.extract_text()
    return text

def evaluate_requirements(policy_text, submission_text):
    prompt = f"""Compare the following two documents and determine if the second document meets the requirements specified in the first document.
    Policy Document: {policy_text}
    Submitted Document: {submission_text}
    Respond with exactly one of these words: GREEN, YELLOW, RED followed by explanation if YELLOW."""
    
    response = bedrock.invoke_model(
        modelId="anthropic.claude-v2",
        body=json.dumps({
            "prompt": prompt,
            "max_tokens_to_sample": 500,
            "temperature": 0
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
