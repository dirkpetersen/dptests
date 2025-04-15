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
from textractor import Textractor
from textractor.data.constants import TextractFeatures
from textractor.data.text_linearization_config import TextLinearizationConfig

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize Textractor
textractor = Textractor(region_name="us-west-2")

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

def process_textract_response(document, config=None):
    """Process Textractor document with optional linearization config"""
    try:
        return document.get_text(config=config) if config else document.get_text()
    except Exception as e:
        logger.error(f"Text processing failed: {str(e)}")
        raise

def process_textract_blocks(blocks):
    """Convert Textract blocks to markdown with tables and layout"""
    output = []
    current_table = []
    in_table = False

    for block in blocks:
        if block['BlockType'] == 'LAYOUT':
            layout_type = block.get('LayoutType', '')
            text = block.get('Text', '').strip()
            
            if layout_type == 'TITLE':
                output.append(f"# {text}\n\n")
            elif layout_type == 'SECTION_HEADER':
                output.append(f"## {text}\n\n")
            elif layout_type == 'HEADER':
                output.append(f"### {text}\n\n")
            elif layout_type == 'FOOTER':
                continue  # Skip footers
            else:
                output.append(f"{text}\n\n")

        elif block['BlockType'] == 'TABLE':
            if not in_table:
                current_table = []
                in_table = True
            # Tables are processed after all their cells are collected
            continue

        elif block['BlockType'] == 'CELL':
            if in_table:
                current_table.append({
                    'text': block.get('Text', '').strip(),
                    'row': block['RowIndex'],
                    'col': block['ColumnIndex'],
                    'header': block.get('Header', False)
                })

        elif block['BlockType'] == 'LINE':
            if in_table:
                # Finish current table
                if current_table:
                    output.append(format_table(current_table))
                current_table = []
                in_table = False
            output.append(f"{block['Text']}\n\n")

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

def textract_async(file_name, bucket_name, region):
    """Extract text from PDF in S3 using async Textract with Textractor"""
    try:
        # Update region in case it's different from initialization
        textractor.region_name = region
        
        document = textractor.analyze_document(
            file_source=f"s3://{bucket_name}/{file_name}",
            features=[TextractFeatures.LAYOUT, TextractFeatures.TABLES],
            save_image=False
        )
        return process_textract_response(document)
    except Exception as e:
        logger.error(f"Async Textract failed: {str(e)}")
        raise

def process_pdf(input_path, output_path, s3_bucket=None, aws_region='us-west-2', max_workers=32, save_images=False, markdown=False):
    """Process a single PDF file using Textract and save as Markdown"""
    config = get_linearization_config(markdown)
    try:
        # First try processing via images
        logger.info("Attempting image-based processing")
        if process_pdf_via_images(input_path, output_path, aws_region=aws_region, max_workers=max_workers, save_images=save_images, markdown=markdown):
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
                    response = textract.analyze_document(
                        Document={'Bytes': document.read()},
                        FeatureTypes=['TABLES', 'LAYOUT']
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

            # Process blocks into markdown
            text = process_textract_blocks(response['Blocks'])

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

def process_image_with_textract(image, region='us-west-2', config=None):
    """Process a single image with Textractor including layout and tables"""
    try:
        # Update region in case it's different from initialization
        textractor.region_name = region
        
        # Convert image to bytes
        img_byte_arr = io.BytesIO()
        image.save(img_byte_arr, format='PNG')
        img_byte_arr = img_byte_arr.getvalue()
        
        document = textractor.analyze_document(
            file_source=img_byte_arr,
            features=[TextractFeatures.LAYOUT, TextractFeatures.TABLES],
            save_image=False
        )
        return process_textract_response(document, config)
    except Exception as e:
        logger.error(f"Image processing failed: {str(e)}")
        raise

def get_linearization_config(markdown=False):
    """Create TextLinearizationConfig from command-line arguments"""
    if markdown:
        return TextLinearizationConfig(
            hide_header_layout=True,
            table_prefix="\n\n",
            table_suffix="\n\n",
            table_row_prefix="| ",
            table_row_suffix=" |",
            table_column_separator="|",
            table_header_divider="---",
            table_use_header_divider=True
        )
    return None

def process_single_image(args):
    """Helper function to process a single image with its index"""
    i, img, aws_region, config = args
    try:
        page_text = process_image_with_textract(img, region=aws_region, config=config)
        return (i, f"--- PAGE {i+1} ---\n\n{page_text}\n\n")
    except Exception as e:
        logger.error(f"Failed processing page {i+1}: {str(e)}")
        return (i, f"--- PAGE {i+1} ---\n\n[ERROR PROCESSING PAGE]\n\n")
    finally:
        del img

def process_pdf_via_images(input_path, output_path, aws_region='us-west-2', dpi=300, max_workers=32, save_images=False, markdown=False):
    """Process PDF by converting to images first with parallel execution"""
    config = get_linearization_config(markdown)
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
            
            # Create list of arguments with indices and config
            tasks = [(i, img, aws_region, config) for i, img in enumerate(images)]
            
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

def convert_folder(input_dir, output_dir, s3_bucket=None, aws_region='us-west-2', max_workers=5, save_images=False, markdown=False):
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
            if process_pdf(input_path, output_path, s3_bucket, aws_region, max_workers, save_images, markdown):
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
    parser.add_argument('--workers', type=int, default=32, 
                      help='Max parallel workers for image processing (1-32, default: 32)')
    parser.add_argument('--images', action='store_true',
                      help='Save page images to output directory')
    parser.add_argument('--markdown', action='store_true',
                      help='Use optimized Markdown formatting')
    
    args = parser.parse_args()
    
    # Validate worker count
    if args.workers < 1 or args.workers > 32:
        logger.warning(f"Invalid worker count {args.workers}, using clamped value")
        args.workers = max(1, min(args.workers, 32))
    
    # Set output directory to input directory if not provided
    if not args.output:
        args.output = args.input
    
    convert_folder(args.input, args.output, args.s3_bucket, args.aws_region, args.workers, args.images, args.markdown)

if __name__ == "__main__":
    main()
