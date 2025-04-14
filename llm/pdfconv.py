import os
import boto3
import argparse
import logging
import PyPDF2
from botocore.exceptions import ClientError

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def process_pdf(input_path, output_path):
    """Process a single PDF file using Textract and save as Markdown"""
    textract = boto3.client('textract')
    
    try:
        # Check for encrypted PDF first
        with open(input_path, 'rb') as f:
            pdf_reader = PyPDF2.PdfReader(f)
            if pdf_reader.is_encrypted:
                logger.error(f"Skipping encrypted PDF: {input_path}")
                return False
                
    except Exception as e:
        logger.error(f"Invalid PDF file {input_path}: {str(e)}")
        return False

    try:
        with open(input_path, 'rb') as document:
            response = textract.detect_document_text(
                Document={'Bytes': document.read()}
            )
    except ClientError as e:
        error_code = e.response['Error']['Code']
        if error_code == 'UnsupportedDocumentException':
            logger.error(f"Unsupported PDF format in {input_path} - file may be corrupted or contain non-standard encoding")
        else:
            logger.error(f"AWS Error processing {input_path}: {e}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error processing {input_path}: {e}")
        return False

    # Extract text from response
    text = ""
    for item in response.get("Blocks", []):
        if item["BlockType"] == "LINE":
            text += item["Text"] + "\n\n"

    # Write to markdown file
    with open(output_path, 'w', encoding='utf-8') as md_file:
        md_file.write(text)
    
    return True

def convert_folder(input_dir, output_dir):
    """Convert all PDFs in a folder to markdown files"""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logger.info(f"Created output directory: {output_dir}")

    converted = 0
    for filename in os.listdir(input_dir):
        if filename.lower().endswith('.pdf'):
            input_path = os.path.join(input_dir, filename)
            base_name = os.path.splitext(filename)[0]
            output_path = os.path.join(output_dir, f"{base_name}.md")
            
            logger.info(f"Converting {filename}...")
            if process_pdf(input_path, output_path):
                converted += 1
            else:
                logger.warning(f"Failed to convert {filename}")

    logger.info(f"Conversion complete. {converted} files converted.")

def main():
    parser = argparse.ArgumentParser(description='Convert PDF files to Markdown using AWS Textract')
    parser.add_argument('-i', '--input', required=True, help='Input directory containing PDF files')
    parser.add_argument('-o', '--output', required=False, help='Output directory for Markdown files (default: same as input directory)')
    
    args = parser.parse_args()
    
    # Set output directory to input directory if not provided
    if not args.output:
        args.output = args.input
    
    convert_folder(args.input, args.output)

if __name__ == "__main__":
    main()
