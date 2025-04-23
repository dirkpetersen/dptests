import os
import boto3
import argparse
import logging
import fitz  # PyMuPDF
import io
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

def process_textract_blocks(blocks):
    """Convert Textract blocks to markdown with tables"""
    output = []
    current_table = []
    in_table = False

    for block in blocks:
        if block['BlockType'] == 'LINE':
            if in_table:
                # Finish current table
                if current_table:
                    output.append(format_table(current_table))
                current_table = []
                in_table = False
            output.append(f"{block['Text']}\n\n")
        
        elif block['BlockType'] == 'TABLE':
            in_table = True
            # Tables are processed after all their cells are collected
            continue
        
        elif block['BlockType'] == 'CELL':
            if in_table:
                current_table.append({
                    'text': block.get('Text', '').strip(),
                    'row': block['RowIndex'],
                    'col': block['ColumnIndex'],
                    'header': 'COLUMN_HEADER' in block.get('EntityTypes', [])
                })

    # Add any remaining table
    if current_table:
        output.append(format_table(current_table))

    return ''.join(output)

def format_table(cells):
    """Convert table cells to markdown table format"""
    # Find table dimensions
    max_row = max(c['row'] for c in cells)
    max_col = max(c['col'] for c in cells)
    
    # Create empty table structure
    table = [[None for _ in range(max_col)] for _ in range(max_row)]
    headers = set()
    
    # Populate table
    for cell in cells:
        row = cell['row'] - 1  # Textract uses 1-based indexing
        col = cell['col'] - 1
        table[row][col] = cell['text']
        if cell['header']:
            headers.add(row)
    
    # Build markdown table
    md_table = []
    for i, row in enumerate(table):
        md_table.append('| ' + ' | '.join(cell if cell else '' for cell in row) + ' |')
        if i == 0 or i in headers:
            md_table.append('|' + '|'.join(['---'] * len(row)) + '|')
    
    return '\n'.join(md_table) + '\n\n'


def process_pdf(input_path, output_path, s3_bucket=None, aws_region='us-west-2', max_workers=32, save_images=False):
    """Process a single PDF file using Textract and save as Markdown"""
    textract = boto3.client('textract', region_name=aws_region)
    
    try:
        # First try direct PDF processing
        with open(input_path, 'rb') as file:
            img_byte_arr = file.read()

        # Check file size (Textract sync limit is 10MB)
        file_size = len(img_byte_arr)
        if file_size > 10 * 1024 * 1024:
            if not s3_bucket:
                return try_fallback_extraction(input_path, output_path)
            
            # Use async processing for large files
            s3_object = upload_to_s3(input_path, s3_bucket)
            if not s3_object:
                return try_fallback_extraction(input_path, output_path)

            try:
                response = textract.start_document_analysis(
                    DocumentLocation={'S3Object': {'Bucket': s3_bucket, 'Name': s3_object}},
                    FeatureTypes=['TABLES', 'FORMS']
                )
                job_id = response['JobId']
                
                while True:
                    status = textract.get_document_analysis(JobId=job_id)
                    if status['JobStatus'] in ['SUCCEEDED', 'FAILED']:
                        break
                    time.sleep(5)
                
                if status['JobStatus'] != 'SUCCEEDED':
                    raise Exception("Textract job failed")
                
                blocks = []
                next_token = None
                while True:
                    response = textract.get_document_analysis(
                        JobId=job_id,
                        NextToken=next_token
                    ) if next_token else textract.get_document_analysis(JobId=job_id)
                    
                    blocks.extend(response['Blocks'])
                    next_token = response.get('NextToken')
                    if not next_token:
                        break
                
                text = process_textract_blocks(blocks)
                with open(output_path, 'w', encoding='utf-8') as f:
                    f.write(text)
                return True
            finally:
                delete_from_s3(s3_bucket, s3_object)
        else:
            # Sync processing for small files
            response = textract.analyze_document(
                Document={'Bytes': img_byte_arr},
                FeatureTypes=['TABLES', 'FORMS']
            )
            text = process_textract_blocks(response['Blocks'])
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(text)
            return True

    except Exception as e:
        logger.error(f"Processing failed: {str(e)}")
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
    """Process a single image with Textract including tables"""
    try:
        textract = boto3.client('textract', region_name=region)
        
        # Convert image to bytes
        img_byte_arr = io.BytesIO()
        image.save(img_byte_arr, format='PNG')
        img_byte_arr = img_byte_arr.getvalue()
        
        response = textract.analyze_document(
            Document={'Bytes': img_byte_arr},
            FeatureTypes=['TABLES', 'FORMS']
        )
        
        return process_textract_blocks(response['Blocks'])
    except Exception as e:
        logger.error(f"Image processing failed: {str(e)}")
        raise

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

def process_pdf_via_images(input_path, output_path, aws_region='us-west-2', dpi=300, max_workers=32, save_images=False):
    """Process PDF by converting to images first with parallel execution"""
    try:
        logger.info('Converting PDF to images...')
        images = convert_pdf_to_images(input_path, dpi)
        full_text = []
        
        # Create base filename from input PDF
        base_name = os.path.splitext(os.path.basename(input_path))[0]
        output_dir = os.path.dirname(output_path)
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Save images if requested
            if save_images:
                for i, img in enumerate(images):
                    img_path = os.path.join(output_dir, f"{base_name}_page{i+1}.png")
                    img.save(img_path, "PNG")
                    logger.info(f"Saved image: {img_path}")
            # Calculate optimal worker count
            num_pages = len(images)
            num_workers = min(num_pages, max_workers)
            num_workers = max(num_workers, 1)  # Ensure at least 1 worker
            
            # Create list of arguments with indices
            tasks = [(i, img, aws_region) for i, img in enumerate(images)]
            
            logger.info(f"Processing {num_pages} pages with {num_workers} workers")
            
            with ThreadPoolExecutor(max_workers=num_workers) as executor:
                futures = [executor.submit(process_single_image, task) for task in tasks]
                
                for future in as_completed(futures):
                    try:
                        page_num, text = future.result()
                        full_text.append((page_num, text))
                        if (page_num + 1) % 5 == 0:
                            logger.info(f"Completed {page_num + 1}/{num_pages} pages")
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

def convert_folder(input_dir, output_dir, s3_bucket=None, aws_region='us-west-2', max_workers=5, save_images=False):
    """Convert all PDFs in a folder to markdown files"""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logger.info(f"Created output directory: {output_dir}")

    converted = 0
    pdf_files = [f for f in os.listdir(input_dir) if f.lower().endswith('.pdf')]
    for filename in pdf_files:
        if filename.lower().endswith('.pdf'):
            input_path = os.path.join(input_dir, filename)
            base_name = os.path.splitext(filename)[0]
            output_path = os.path.join(output_dir, f"{base_name}.md")
            
            logger.info(f"Converting {filename}...")
            try:
                if save_images:
                    # Use image-based processing path
                    success = process_pdf_via_images(
                        input_path, 
                        output_path,
                        aws_region=aws_region,
                        dpi=300,
                        max_workers=max_workers,
                        save_images=True
                    )
                else:
                    # Use direct PDF processing
                    success = process_pdf(
                        input_path,
                        output_path,
                        s3_bucket,
                        aws_region,
                        max_workers,
                        save_images=False
                    )
                
                if success:
                    converted += 1
                else:
                    logger.warning(f"Failed to convert {filename}")
            except Exception as e:
                logger.error(f"Error processing {filename}: {str(e)}")
                continue

    logger.info(f"Conversion complete. {converted}/{len(pdf_files)} files converted.")

def main():
    parser = argparse.ArgumentParser(description='Convert PDF files to Markdown using AWS Textract')
    parser.add_argument('-i', '--input', required=True, help='Input directory containing PDF files')
    parser.add_argument('-o', '--output', required=False, help='Output directory for Markdown files (default: same as input directory)')
    parser.add_argument('--s3-bucket', help='S3 bucket for large files (>10MB)')
    parser.add_argument('--aws-region', default='us-west-2', help='AWS region')
    parser.add_argument('--workers', type=int, default=32, 
                      help='Max parallel workers for image processing (1-32, default: 32)')
    parser.add_argument('--images', action='store_true',
                      help='Save page images to output directory')
    
    args = parser.parse_args()
    
    # Validate worker count
    if args.workers < 1 or args.workers > 32:
        logger.warning(f"Invalid worker count {args.workers}, using clamped value")
        args.workers = max(1, min(args.workers, 32))
    
    # Set output directory to input directory if not provided
    if not args.output:
        args.output = args.input
    
    convert_folder(args.input, args.output, args.s3_bucket, args.aws_region, args.workers, args.images)

if __name__ == "__main__":
    main()
