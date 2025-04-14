import os
import boto3
import argparse
import logging
import PyPDF2
import fitz  # PyMuPDF
import io
import shutil
import time
from PIL import Image
import tempfile
from botocore.exceptions import ClientError
from concurrent.futures import ThreadPoolExecutor, as_completed

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def upload_to_s3(file_path, bucket_name, object_name=None):
    """Upload a file to an S3 bucket"""
    s3 = boto3.client('s3')
    if object_name is None:
        object_name = os.path.basename(file_path)
    try:
        s3.upload_file(file_path, bucket_name, object_name)
        return object_name
    except ClientError as e:
        logger.error(f"Failed to upload {file_path} to S3: {e}")
        return None

def delete_from_s3(bucket_name, object_name):
    """Delete an object from S3"""
    s3 = boto3.client('s3')
    try:
        s3.delete_object(Bucket=bucket_name, Key=object_name)
    except ClientError as e:
        logger.error(f"Failed to delete {object_name} from S3: {e}")

def textract_async(file_name, bucket_name, region):
    """Extract text from PDF in S3 using async Textract"""
    textract_client = boto3.client('textract', region_name=region)
    
    response = textract_client.start_document_text_detection(
        DocumentLocation={'S3Object': {'Bucket': bucket_name, 'Name': file_name}}
    )
    
    job_id = response['JobId']
    logger.info(f"Started Textract job {job_id}")
    
    while True:
        response = textract_client.get_document_text_detection(JobId=job_id)
        status = response['JobStatus']
        if status in ['SUCCEEDED', 'FAILED']:
            break
        time.sleep(5)
    
    if status != 'SUCCEEDED':
        raise Exception(f"Textract job failed: {status}")
    
    pages = []
    next_token = None
    
    while True:
        if next_token:
            response = textract_client.get_document_text_detection(
                JobId=job_id,
                NextToken=next_token
            )
        else:
            response = textract_client.get_document_text_detection(JobId=job_id)
        
        pages.extend(response['Blocks'])
        next_token = response.get('NextToken')
        if not next_token:
            break
    
    text = ""
    for item in pages:
        if item["BlockType"] == "LINE":
            text += item["Text"] + "\n\n"
    
    return text

def process_pdf(input_path, output_path, s3_bucket=None, aws_region='us-west-2', max_workers=5):
    """Process a single PDF file using Textract and save as Markdown"""
    try:
        # First try processing via images
        logger.info("Attempting image-based processing")
        if process_pdf_via_images(input_path, output_path, aws_region=aws_region, max_workers=max_workers):
            return True
        
        # If image processing failed, fall back to previous methods
        logger.info("Falling back to direct PDF processing")
        
        # Check for encrypted PDF first
        with open(input_path, 'rb') as f:
            pdf_reader = PyPDF2.PdfReader(f)
            if pdf_reader.is_encrypted:
                logger.error(f"Skipping encrypted PDF: {input_path}")
                return False
        
        # Check file size (Textract sync limit is 10MB)
        file_size = os.path.getsize(input_path)
        if file_size > 10 * 1024 * 1024:
            if not s3_bucket:
                logger.warning(f"File too large for sync processing (10MB max) and no S3 bucket provided: {input_path}")
                return try_fallback_extraction(input_path, output_path)
            
            # Process large files with async Textract via S3
            logger.info(f"Processing large file with async Textract: {input_path}")
            
            # Preprocess PDF with PyMuPDF
            doc = fitz.open(input_path)
            temp_path = os.path.join("/tmp", os.path.basename(input_path))
            doc.save(temp_path, deflate=True, garbage=4)
            doc.close()
            
            # Upload to S3
            s3_object = upload_to_s3(temp_path, s3_bucket)
            if not s3_object:
                os.remove(temp_path)
                return try_fallback_extraction(input_path, output_path)
                
            try:
                # Process with async Textract
                text = textract_async(s3_object, s3_bucket, aws_region)
                
                # Write to markdown file
                with open(output_path, 'w', encoding='utf-8') as md_file:
                    md_file.write(text)
                
                return True
            finally:
                # Clean up
                delete_from_s3(s3_bucket, s3_object)
                if os.path.exists(temp_path):
                    os.remove(temp_path)
        else:
            # For small files, use sync Textract API
            textract = boto3.client('textract')
            
            # Preprocess PDF with PyMuPDF
            doc = fitz.open(input_path)
            pdf_bytes = io.BytesIO()
            doc.save(pdf_bytes, deflate=True, garbage=4)
            pdf_bytes.seek(0)
            doc.close()

            # Create temp processed PDF
            temp_path = os.path.join("/tmp", os.path.basename(input_path))
            with open(temp_path, "wb") as f:
                f.write(pdf_bytes.getbuffer())
                    
            try:
                with open(temp_path, 'rb') as document:
                    response = textract.detect_document_text(
                        Document={'Bytes': document.read()}
                    )
            except ClientError as e:
                error_code = e.response['Error']['Code']
                if error_code == 'UnsupportedDocumentException':
                    logger.error(f"Unsupported PDF format in {input_path} - trying fallback extraction")
                    return try_fallback_extraction(input_path, output_path)
                else:
                    logger.error(f"AWS Error processing {input_path}: {e}")
                    return False
            except Exception as e:
                logger.error(f"Unexpected error processing {input_path}: {e}")
                return False
            finally:
                if os.path.exists(temp_path):
                    os.remove(temp_path)

            # Extract text from response
            text = ""
            for item in response.get("Blocks", []):
                if item["BlockType"] == "LINE":
                    text += item["Text"] + "\n\n"

            # Write to markdown file
            with open(output_path, 'w', encoding='utf-8') as md_file:
                md_file.write(text)
            
            return True
                
    except Exception as e:
        logger.error(f"Invalid PDF file {input_path}: {str(e)}")
        return try_fallback_extraction(input_path, output_path)

def convert_pdf_to_images(pdf_path, dpi=300):
    """Convert PDF to list of PIL Images"""
    doc = fitz.open(pdf_path)
    images = []
    
    for page in doc:
        pix = page.get_pixmap(matrix=fitz.Matrix(dpi/72, dpi/72))
        img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
        images.append(img)
    
    return images

def process_image_with_textract(image, region='us-west-2'):
    """Process a single image with Textract"""
    textract = boto3.client('textract', region_name=region)
    
    # Convert image to bytes
    img_byte_arr = io.BytesIO()
    image.save(img_byte_arr, format='PNG')
    img_byte_arr = img_byte_arr.getvalue()
    
    response = textract.detect_document_text(
        Document={'Bytes': img_byte_arr}
    )
    
    text = ""
    for item in response.get("Blocks", []):
        if item["BlockType"] == "LINE":
            text += item["Text"] + "\n\n"
    
    return text

def process_single_image(args):
    """Helper function to process a single image with its index"""
    i, img, aws_region = args
    try:
        page_text = process_image_with_textract(img, region=aws_region)
        return (i, f"--- PAGE {i+1} ---\n\n{page_text}\n\n")
    except Exception as e:
        logger.error(f"Failed processing page {i+1}: {str(e)}")
        return (i, f"--- PAGE {i+1} ---\n\n[ERROR PROCESSING PAGE]\n\n")
    finally:
        del img

def process_pdf_via_images(input_path, output_path, aws_region='us-west-2', dpi=300, max_workers=5):
    """Process PDF by converting to images first with parallel execution"""
    try:
        logger.info('Converting PDF to images...')
        images = convert_pdf_to_images(input_path, dpi)
        full_text = []
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create list of arguments with indices
            tasks = [(i, img, aws_region) for i, img in enumerate(images)]
            
            logger.info(f"Processing {len(images)} pages with {max_workers} workers")
            
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [executor.submit(process_single_image, task) for task in tasks]
                
                for future in as_completed(futures):
                    try:
                        page_num, text = future.result()
                        full_text.append((page_num, text))
                        if (page_num + 1) % 5 == 0:
                            logger.info(f"Completed {page_num + 1}/{len(images)} pages")
                    except Exception as e:
                        logger.error(f"Page processing failed: {str(e)}")

            # Sort results by page number and combine
            full_text.sort(key=lambda x: x[0])
            combined_text = "".join([t[1] for t in full_text])
            
            # Save final text
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(combined_text)
        
        return True
        
    except Exception as e:
        logger.error(f"Image-based processing failed: {str(e)}")
        return False

def try_fallback_extraction(input_path, output_path):
    """Fallback text extraction using PyMuPDF when Textract fails"""
    try:
        doc = fitz.open(input_path)
        text = ""
        for page in doc:
            text += page.get_text("text") + "\n\n"
        
        with open(output_path, 'w', encoding='utf-8') as md_file:
            md_file.write(text)
        return True
    except Exception as e:
        logger.error(f"Fallback extraction failed for {input_path}: {e}")
        return False

def convert_folder(input_dir, output_dir, s3_bucket=None, aws_region='us-west-2', max_workers=5):
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
            if process_pdf(input_path, output_path, s3_bucket, aws_region, max_workers):
                converted += 1
            else:
                logger.warning(f"Failed to convert {filename}")

    logger.info(f"Conversion complete. {converted} files converted.")

def main():
    parser = argparse.ArgumentParser(description='Convert PDF files to Markdown using AWS Textract')
    parser.add_argument('-i', '--input', required=True, help='Input directory containing PDF files')
    parser.add_argument('-o', '--output', required=False, help='Output directory for Markdown files (default: same as input directory)')
    parser.add_argument('--s3-bucket', help='S3 bucket for large files (>10MB)')
    parser.add_argument('--aws-region', default='us-west-2', help='AWS region')
    parser.add_argument('--workers', type=int, default=5, 
                      help='Max parallel workers for image processing (default: 5)')
    
    args = parser.parse_args()
    
    # Set output directory to input directory if not provided
    if not args.output:
        args.output = args.input
    
    convert_folder(args.input, args.output, args.s3_bucket, args.aws_region, args.workers)

if __name__ == "__main__":
    main()
