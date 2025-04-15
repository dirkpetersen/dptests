import argparse
import logging
import pymupdf4llm
from pathlib import Path

def convert_pdf_folder(input_folder: str, output_dir: str = "markdown_output"):
    """Convert all PDFs in a folder to markdown files preserving structure"""
    input_path = Path(input_folder)
    output_path = Path(output_dir)
    
    # Create output directory if needed
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find all PDF files recursively
    pdf_files = list(input_path.glob("**/*.pdf"))
    
    if not pdf_files:
        logging.warning(f"No PDF files found in {input_folder}")
        return

    logging.info(f"Found {len(pdf_files)} PDF files to convert")
    
    for pdf_file in pdf_files:
        try:
            # Generate output path with same relative structure
            relative_path = pdf_file.relative_to(input_path)
            md_path = output_path / relative_path.with_suffix(".md")
            
            # Create parent directories if needed
            md_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Convert and save
            logging.info(f"Converting {pdf_file}...")
            md_text = pymupdf4llm.to_markdown(str(pdf_file))
            md_path.write_text(md_text, encoding="utf-8")
            
        except Exception as e:
            logging.error(f"Failed to convert {pdf_file}: {str(e)}")

def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    
    parser = argparse.ArgumentParser(description="Convert PDF files in a folder to markdown")
    parser.add_argument("input_folder", help="Path to folder containing PDF files")
    parser.add_argument("-o", "--output", default="markdown_output", 
                      help="Output directory for markdown files")
    args = parser.parse_args()
    
    convert_pdf_folder(args.input_folder, args.output)
    logging.info("Conversion complete")

if __name__ == "__main__":
    main()
