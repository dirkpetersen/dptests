#!/usr/bin/python3

import os
import sys
import argparse
from pathlib import Path
import fitz  # PyMuPDF
from tqdm import tqdm
import signal
from contextlib import contextmanager
try:
    import pymupdf4llm
except ImportError:
    pymupdf4llm = None


@contextmanager
def timeout_handler(seconds):
    """
    Context manager to handle timeouts for problematic operations.
    """
    def signal_handler(signum, frame):
        raise TimeoutError("Operation timed out")
    
    # Set up signal handler
    old_handler = signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)


def detect_table_pages(pdf_path):
    """
    Detect which pages contain tables using PyMuPDF's find_tables().
    
    Args:
        pdf_path (str): Path to the PDF file
    
    Returns:
        list: List of page numbers (1-indexed) that contain tables
    """
    try:
        doc = fitz.open(str(pdf_path))
        table_pages = []
        
        try:
            # Check each page for tables with timeout
            for page_num in range(len(doc)):
                try:
                    with timeout_handler(5):  # 5 second timeout per page
                        page = doc[page_num]
                        tables = page.find_tables()
                        
                        if tables:
                            # Verify at least one table has content
                            for table in tables:
                                try:
                                    table_data = table.extract()
                                    if table_data and len(table_data) > 0 and len(table_data[0]) > 0:
                                        table_pages.append(page_num + 1)  # Convert to 1-indexed
                                        break
                                except Exception:
                                    continue
                except (KeyboardInterrupt, SystemExit):
                    raise
                except (TimeoutError, Exception):
                    # Skip pages that cause errors or timeout
                    continue
        finally:
            doc.close()
        
        return sorted(list(set(table_pages)))
    
    except Exception:
        return []


def pdf_to_markdown_simple(pdf_path, output_path=None):
    """
    Convert PDF to markdown using pymupdf4llm (fast, simple conversion).
    
    Args:
        pdf_path (str): Path to the input PDF file
        output_path (str, optional): Path for the output MD file
    
    Returns:
        str: Path to the created markdown file
    """
    if not pymupdf4llm:
        raise ImportError("pymupdf4llm not available, falling back to standard conversion")
    
    pdf_path = Path(pdf_path)
    
    if output_path is None:
        output_path = pdf_path.with_suffix('.md')
    else:
        output_path = Path(output_path)
    
    if not pdf_path.exists():
        raise FileNotFoundError(f"PDF file not found: {pdf_path}")
    
    # Use pymupdf4llm for fast conversion
    md_text = pymupdf4llm.to_markdown(str(pdf_path))
    
    # Write to output file
    output_path.write_bytes(md_text.encode('utf-8'))
    
    return str(output_path)


def pdf_to_markdown(pdf_path, output_path=None, show_progress=True):
    """
    Convert a PDF file to Markdown format with special attention to tables.
    
    Args:
        pdf_path (str): Path to the input PDF file
        output_path (str, optional): Path for the output MD file. If None, uses same base name as PDF.
        show_progress (bool): Whether to show progress bar for page processing
    
    Returns:
        str: Path to the created markdown file
    """
    pdf_path = Path(pdf_path)
    
    if output_path is None:
        output_path = pdf_path.with_suffix('.md')
    else:
        output_path = Path(output_path)
    
    if not pdf_path.exists():
        raise FileNotFoundError(f"PDF file not found: {pdf_path}")
    
    # Table detection is now handled in main() to avoid duplication
    
    doc = fitz.open(str(pdf_path))
    markdown_content = []
    
    try:
        page_range = range(len(doc))
        if show_progress and len(doc) > 1:
            page_range = tqdm(page_range, desc=f"Processing {pdf_path.name}", unit="pages")
        
        for page_num in page_range:
            try:
                page = doc[page_num]
                
                # Extract text with layout information
                blocks = page.get_text("dict")
                
                # Process blocks to identify tables and regular text
                page_content = []
                
                for block in blocks["blocks"]:
                    if "lines" in block:  # Text block
                        block_text = []
                        for line in block["lines"]:
                            line_text = ""
                            for span in line["spans"]:
                                text = span["text"].strip()
                                if text:
                                    line_text += text + " "
                            if line_text.strip():
                                block_text.append(line_text.strip())
                        
                        if block_text:
                            # Check if this looks like a table (multiple aligned columns)
                            if _is_table_block(block_text):
                                table_md = _format_as_table(block_text)
                                page_content.append(table_md)
                            else:
                                page_content.extend(block_text)
                
                # Try to extract actual tables using PyMuPDF's table detection
                try:
                    # Only attempt table detection if the page is not too complex
                    tables = page.find_tables()
                    for table in tables:
                        try:
                            table_data = table.extract()
                            if table_data:
                                table_md = _format_table_data(table_data)
                                page_content.append(table_md)
                        except (KeyboardInterrupt, SystemExit):
                            raise
                        except Exception:
                            # Skip problematic tables but continue processing
                            continue
                except (KeyboardInterrupt, SystemExit):
                    raise
                except Exception:
                    pass  # Fallback to text extraction if table detection fails
            
                if page_content:
                    if page_num > 0:
                        markdown_content.append(f"\n---\n# Page {page_num + 1}\n")
                    markdown_content.extend(page_content)
                    markdown_content.append("")
            
            except (KeyboardInterrupt, SystemExit):
                raise
            except Exception as e:
                # Skip problematic pages but continue processing
                if show_progress:
                    tqdm.write(f"Warning: Skipped page {page_num + 1} due to error: {str(e)[:50]}...")
                continue
    
    finally:
        doc.close()
    
    # Write to output file
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(markdown_content))
    
    return str(output_path)


def _is_table_block(text_lines):
    """Check if a block of text looks like a table based on alignment patterns."""
    if len(text_lines) < 2:
        return False
    
    # Look for consistent spacing patterns that suggest columns
    space_patterns = []
    for line in text_lines[:5]:  # Check first 5 lines
        spaces = []
        for i, char in enumerate(line):
            if char == ' ' and (i == 0 or line[i-1] != ' '):
                spaces.append(i)
        space_patterns.append(spaces)
    
    # If multiple lines have similar space patterns, likely a table
    if len(space_patterns) >= 2:
        first_pattern = space_patterns[0]
        similar_count = sum(1 for pattern in space_patterns[1:] 
                          if len(set(first_pattern) & set(pattern)) >= len(first_pattern) // 2)
        return similar_count >= len(space_patterns) // 2
    
    return False


def _format_as_table(text_lines):
    """Format aligned text as a markdown table."""
    if not text_lines:
        return ""
    
    # Split lines into columns based on multiple spaces
    rows = []
    for line in text_lines:
        # Split on 2+ spaces to identify columns
        columns = [col.strip() for col in line.split('  ') if col.strip()]
        if columns:
            rows.append(columns)
    
    if not rows:
        return '\n'.join(text_lines)
    
    # Determine max columns
    max_cols = max(len(row) for row in rows)
    
    # Normalize all rows to have the same number of columns
    normalized_rows = []
    for row in rows:
        normalized_row = row + [''] * (max_cols - len(row))
        normalized_rows.append(normalized_row)
    
    # Create markdown table
    if normalized_rows:
        table_md = []
        # Header row
        table_md.append('| ' + ' | '.join(normalized_rows[0]) + ' |')
        # Separator row
        table_md.append('| ' + ' | '.join(['---'] * max_cols) + ' |')
        # Data rows
        for row in normalized_rows[1:]:
            table_md.append('| ' + ' | '.join(row) + ' |')
        
        return '\n' + '\n'.join(table_md) + '\n'
    
    return '\n'.join(text_lines)


def _format_table_data(table_data):
    """Format extracted table data as markdown table."""
    if not table_data or not table_data[0]:
        return ""
    
    table_md = []
    
    # Header row
    header = [str(cell) if cell else '' for cell in table_data[0]]
    table_md.append('| ' + ' | '.join(header) + ' |')
    
    # Separator row
    table_md.append('| ' + ' | '.join(['---'] * len(header)) + ' |')
    
    # Data rows
    for row in table_data[1:]:
        formatted_row = [str(cell) if cell else '' for cell in row]
        # Pad row to match header length
        while len(formatted_row) < len(header):
            formatted_row.append('')
        table_md.append('| ' + ' | '.join(formatted_row) + ' |')
    
    return '\n' + '\n'.join(table_md) + '\n'


def merge_markdown_files(md_files, output_path, folder_name):
    """
    Merge multiple markdown files into a single document.
    
    Args:
        md_files (list): List of markdown file paths
        output_path (str): Path for the merged output file
        folder_name (str): Name of the source folder for the title
    
    Returns:
        str: Path to the merged markdown file
    """
    merged_content = []
    
    for md_file in md_files:
        md_path = Path(md_file)
        original_pdf_name = md_path.stem + '.pdf'
        
        # Add filename header
        merged_content.append(f"# {original_pdf_name}\n")
        
        # Read and add file content
        try:
            with open(md_path, 'r', encoding='utf-8') as f:
                content = f.read().strip()
                if content:
                    merged_content.append(content)
                    merged_content.append("\n\n---\n")  # Separator between files
        except Exception as e:
            merged_content.append(f"Error reading {original_pdf_name}: {e}\n\n---\n")
    
    # Remove the last separator
    if merged_content and merged_content[-1] == "\n\n---\n":
        merged_content.pop()
    
    # Write merged content
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(merged_content))
    
    return str(output_path)


def process_folder(folder_path, output_folder=None, show_progress=True, show_table_detection=True, use_fast_conversion=False, merge=False):
    """
    Convert all PDF files in a folder to markdown.
    
    Args:
        folder_path (str): Path to folder containing PDF files
        output_folder (str, optional): Output folder path. If None, uses same folder.
        show_progress (bool): Whether to show progress bar for file processing
        show_table_detection (bool): Whether to show table detection for each file
        use_fast_conversion (bool): Whether to use pymupdf4llm for fast conversion
        merge (bool): Whether to merge all markdown files into a single document
    
    Returns:
        list: Paths to created markdown files (or single merged file if merge=True)
    """
    folder_path = Path(folder_path)
    
    if not folder_path.exists() or not folder_path.is_dir():
        raise ValueError(f"Folder not found or not a directory: {folder_path}")
    
    if output_folder:
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)
    else:
        output_folder = folder_path
    
    pdf_files = list(folder_path.glob("*.pdf"))
    
    if not pdf_files:
        print(f"No PDF files found in {folder_path}")
        return []
    
    converted_files = []
    
    file_iterator = pdf_files
    if show_progress:
        file_iterator = tqdm(pdf_files, desc="Converting PDFs", unit="files")
    
    for pdf_file in file_iterator:
        try:
            # Show table detection for each file if enabled
            if show_table_detection:
                table_pages = detect_table_pages(pdf_file)
                if table_pages:
                    if show_progress:
                        tqdm.write(f"{pdf_file.name}: Found tables in pages {', '.join(map(str, table_pages))}")
                    else:
                        print(f"{pdf_file.name}: Found tables in pages {', '.join(map(str, table_pages))}")
            
            output_path = output_folder / (pdf_file.stem + '.md')
            
            # Use fast conversion if enabled
            if use_fast_conversion:
                try:
                    result_path = pdf_to_markdown_simple(pdf_file, output_path)
                except Exception:
                    # Fallback to standard conversion
                    result_path = pdf_to_markdown(pdf_file, output_path, show_progress=False)
            else:
                result_path = pdf_to_markdown(pdf_file, output_path, show_progress=False)
            
            converted_files.append(result_path)
            if not show_progress:
                print(f"Converted: {pdf_file.name} -> {Path(result_path).name}")
        except Exception as e:
            if show_progress:
                tqdm.write(f"Error converting {pdf_file.name}: {e}")
            else:
                print(f"Error converting {pdf_file.name}: {e}")
    
    # If merge is enabled, merge all markdown files into one
    if merge and converted_files:
        merged_filename = f"{folder_path.name}.md"
        merged_path = output_folder / merged_filename
        
        # Sort files by name for consistent ordering
        sorted_files = sorted(converted_files, key=lambda x: Path(x).name.lower())
        
        merged_file = merge_markdown_files(sorted_files, merged_path, folder_path.name)
        
        # Clean up individual markdown files
        for md_file in converted_files:
            try:
                Path(md_file).unlink()
            except Exception:
                pass  # Continue if cleanup fails
        
        print(f"Merged {len(converted_files)} files into: {merged_filename}")
        return [merged_file]
    
    return converted_files


def main():
    parser = argparse.ArgumentParser(
        description="Convert PDF files to Markdown format with table extraction",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s document.pdf                    # Convert to document.md
  %(prog)s document.pdf output.md          # Convert to output.md
  %(prog)s /path/to/pdfs/                  # Convert all PDFs in folder
  %(prog)s /path/to/pdfs/ /path/to/output/ # Convert all PDFs to output folder
  %(prog)s /path/to/pdfs/ --merge          # Merge all PDFs into single [directory-name].md
        """
    )
    
    parser.add_argument('input', help='PDF file or folder containing PDF files')
    parser.add_argument('output', nargs='?', help='Output file/folder (optional)')
    parser.add_argument('--no-progress', action='store_true', help='Disable progress bars')
    parser.add_argument('--no-table-detection', action='store_true', help='Skip table detection phase')
    parser.add_argument('--merge', action='store_true', help='Merge all PDFs in folder into single markdown file named [directory-name].md')
    
    args = parser.parse_args()
    
    input_path = Path(args.input)
    
    try:
        if input_path.is_file():
            # Single file conversion
            if input_path.suffix.lower() != '.pdf':
                print(f"Error: Input file must be a PDF file, got: {input_path.suffix}")
                sys.exit(1)
            
            if args.merge:
                print("Error: --merge flag can only be used with folders, not individual files")
                sys.exit(1)
            
            show_progress = not args.no_progress
            show_table_detection = not args.no_table_detection
            
            # Use fast pymupdf4llm conversion if both progress and table detection are disabled
            if not show_progress and not show_table_detection and pymupdf4llm:
                try:
                    output_path = pdf_to_markdown_simple(input_path, args.output)
                    print(f"Successfully converted: {input_path.name} -> {Path(output_path).name}")
                except Exception as e:
                    print(f"Fast conversion failed, falling back to standard method: {e}")
                    # Fallback to standard conversion
                    if show_table_detection:
                        table_pages = detect_table_pages(input_path)
                        if table_pages:
                            print(f"Found tables in pages: {', '.join(map(str, table_pages))}")
                    output_path = pdf_to_markdown(input_path, args.output, show_progress)
                    print(f"Successfully converted: {input_path.name} -> {Path(output_path).name}")
            else:
                # Standard conversion with table detection and/or progress
                if show_table_detection:
                    table_pages = detect_table_pages(input_path)
                    if table_pages:
                        print(f"Found tables in pages: {', '.join(map(str, table_pages))}")
                output_path = pdf_to_markdown(input_path, args.output, show_progress)
                print(f"Successfully converted: {input_path.name} -> {Path(output_path).name}")
        
        elif input_path.is_dir():
            # Folder conversion
            show_progress = not args.no_progress
            show_table_detection = not args.no_table_detection
            use_fast_conversion = not show_progress and not show_table_detection and pymupdf4llm
            converted_files = process_folder(input_path, args.output, show_progress, show_table_detection, use_fast_conversion, args.merge)
            
            if args.merge:
                if converted_files:
                    # Count original PDFs by checking how many were processed
                    pdf_count = len(list(input_path.glob("*.pdf")))
                    print(f"Successfully merged {pdf_count} PDF files into: {Path(converted_files[0]).name}")
                else:
                    print("No PDF files found to merge")
            else:
                print(f"Converted {len(converted_files)} PDF files to Markdown")
        
        else:
            print(f"Error: Input path does not exist: {input_path}")
            sys.exit(1)
    
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()