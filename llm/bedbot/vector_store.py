"""
Qdrant Vector Store for BedBot - Document Embeddings and Retrieval with Unstructured.io
"""

import os
import json
import logging
import tempfile
import hashlib
import io
from datetime import datetime
from typing import List, Dict, Tuple, Optional, Any
from pathlib import Path

# Set up logger
logger = logging.getLogger(__name__)

try:
    # Import Qdrant client
    from qdrant_client import QdrantClient
    from qdrant_client.http import models
    from qdrant_client.http.models import Distance, VectorParams, PointStruct
    
    # Import embeddings
    from sentence_transformers import SentenceTransformer
    
    # Import unstructured.io for document processing
    from unstructured.partition.pdf import partition_pdf
    from unstructured.staging.base import elements_to_json
    from unstructured.chunking.title import chunk_by_title
    
    VECTOR_STORE_AVAILABLE = True  # Renamed from VECTOR_STORE_AVAILABLE for generic compatibility
    logger.info("Qdrant, SentenceTransformers, and Unstructured loaded successfully")
except ImportError as e:
    VECTOR_STORE_AVAILABLE = False  # Renamed from VECTOR_STORE_AVAILABLE for generic compatibility
    logger.warning(f"Qdrant/Unstructured dependencies not available: {e}")
    logger.warning("To enable vector store: pip install qdrant-client sentence-transformers unstructured[pdf]")

class DocumentChunk:
    """Represents a chunk of text from a document"""
    def __init__(self, text: str, source: str, chunk_id: str, metadata: Dict = None):
        self.text = text
        self.source = source  # filename
        self.chunk_id = chunk_id
        self.metadata = metadata or {}
        self.timestamp = datetime.now().isoformat()

class QdrantVectorStore:
    """Qdrant-based vector store for document embeddings and retrieval with unstructured.io"""
    
    def __init__(self, session_id: str, embedding_model: str = "all-MiniLM-L6-v2"):
        """
        Initialize Qdrant vector store
        
        Args:
            session_id: Unique session identifier
            embedding_model: SentenceTransformer model name
        """
        if not VECTOR_STORE_AVAILABLE:
            raise ImportError("Qdrant/Unstructured dependencies not available. Install with: pip install qdrant-client sentence-transformers unstructured[pdf]")
        
        self.session_id = session_id
        self.embedding_model_name = embedding_model
        self.collection_name = f"bedbot_session_{session_id.replace('-', '_')}"
        
        # Initialize embedding model
        try:
            self.embedding_model = SentenceTransformer(embedding_model)
            self.embedding_dim = self.embedding_model.get_sentence_embedding_dimension()
            logger.info(f"Loaded embedding model: {embedding_model} (dim: {self.embedding_dim})")
        except Exception as e:
            logger.error(f"Failed to load embedding model {embedding_model}: {e}")
            raise
        
        # Initialize Qdrant client (in-memory)
        try:
            self.client = QdrantClient(":memory:")
            logger.info(f"Initialized in-memory Qdrant client for session: {session_id}")
        except Exception as e:
            logger.error(f"Failed to initialize Qdrant client: {e}")
            raise
        
        # Create collection
        try:
            self.client.create_collection(
                collection_name=self.collection_name,
                vectors_config=VectorParams(
                    size=self.embedding_dim,
                    distance=Distance.COSINE
                )
            )
            logger.info(f"Created Qdrant collection: {self.collection_name}")
        except Exception as e:
            logger.error(f"Failed to create Qdrant collection: {e}")
            raise
        
        # Document tracking
        self.doc_metadata: Dict[str, Any] = {}
        self.documents: List[DocumentChunk] = []
        self.next_point_id = 0
        
        logger.info(f"Initialized Qdrant vector store for session: {session_id}")
    
    def _process_pdf_with_unstructured(self, pdf_input, is_bytes: bool = False) -> str:
        """
        Process PDF using unstructured.io for complete document extraction
        
        Args:
            pdf_input: Path to PDF file or PDF bytes
            is_bytes: True if pdf_input is bytes, False if it's a file path
            
        Returns:
            Processed text content
        """
        try:
            logger.info("Processing PDF with unstructured.io for complete extraction")
            
            if is_bytes:
                # Handle bytes input by saving to temporary file
                with tempfile.NamedTemporaryFile(delete=False, suffix='.pdf') as temp_file:
                    temp_file.write(pdf_input)
                    temp_path = temp_file.name
                
                try:
                    # Process the temporary file with multiple strategies for completeness
                    elements = self._extract_pdf_elements_comprehensive(temp_path)
                finally:
                    # Clean up temporary file
                    os.unlink(temp_path)
            else:
                # Handle file path input
                elements = self._extract_pdf_elements_comprehensive(pdf_input)
            
            logger.info(f"Extracted {len(elements)} elements from PDF")
            
            # Convert elements to comprehensive text with page tracking
            full_text = self._elements_to_comprehensive_text(elements)
            
            logger.info(f"Processed PDF content: {len(full_text)} characters from {len(elements)} elements")
            
            # Validate completeness
            if len(full_text) < 10000:  # NIST docs should be much larger
                logger.warning(f"PDF content seems incomplete ({len(full_text)} chars), trying fallback")
                return self._fallback_pdf_processing(pdf_input, is_bytes)
            
            return full_text
            
        except Exception as e:
            logger.error(f"Error processing PDF with unstructured.io: {e}")
            # Fallback to simple text extraction if unstructured fails
            return self._fallback_pdf_processing(pdf_input, is_bytes)
    
    def _extract_pdf_elements_comprehensive(self, pdf_path: str):
        """Extract PDF elements using multiple strategies for maximum completeness"""
        elements = []
        strategies = [
            # Try hi_res first for best quality
            {"strategy": "hi_res", "infer_table_structure": True},
            # Fallback to auto if hi_res fails
            {"strategy": "auto", "infer_table_structure": True},
            # Fast fallback if others fail
            {"strategy": "fast", "infer_table_structure": False}
        ]
        
        for i, strategy_config in enumerate(strategies):
            try:
                logger.info(f"Trying PDF extraction strategy {i+1}: {strategy_config['strategy']}")
                
                elements = partition_pdf(
                    filename=pdf_path,
                    strategy=strategy_config["strategy"],
                    infer_table_structure=strategy_config["infer_table_structure"],
                    extract_images_in_pdf=False,
                    include_page_breaks=True,  # Ensure page boundaries are preserved
                    multipage_sections=True,   # Handle multi-page sections
                    chunking_strategy=None     # Don't chunk at this level
                )
                
                # Check if we got reasonable content
                total_text = sum(len(str(elem)) for elem in elements)
                logger.info(f"Strategy {strategy_config['strategy']} extracted {len(elements)} elements, {total_text} total chars")
                
                if elements and total_text > 1000:  # Reasonable amount of content
                    logger.info(f"Successfully extracted PDF with strategy: {strategy_config['strategy']}")
                    break
                    
            except Exception as e:
                logger.warning(f"PDF extraction strategy {strategy_config['strategy']} failed: {e}")
                if i == len(strategies) - 1:  # Last strategy failed
                    raise e
                continue
        
        return elements
    
    def _elements_to_comprehensive_text(self, elements) -> str:
        """Convert elements to comprehensive text preserving document structure and completeness"""
        text_content = []
        current_page = None
        page_content = []
        
        # Group content by page to ensure completeness
        page_elements = {}
        
        for element in elements:
            # Get page number
            page_num = getattr(element, 'metadata', {}).get('page_number', 1)
            if page_num not in page_elements:
                page_elements[page_num] = []
            page_elements[page_num].append(element)
        
        # Process pages in order
        for page_num in sorted(page_elements.keys()):
            page_text = f"\n\n=== PAGE {page_num} ===\n\n"
            page_content = []
            
            for element in page_elements[page_num]:
                element_text = str(element).strip()
                if not element_text:
                    continue
                
                # Get element type
                element_type = getattr(element, 'category', 'Text')
                
                # Handle different element types with better preservation
                if element_type in ['Title', 'Header']:
                    # Preserve hierarchical structure
                    page_content.append(f"\n## {element_text}\n")
                elif element_type == 'Table':
                    # Preserve table structure
                    page_content.append(f"\n**TABLE:**\n{element_text}\n")
                elif element_type in ['Text', 'NarrativeText']:
                    # Regular content
                    page_content.append(element_text)
                elif element_type == 'ListItem':
                    # Preserve list structure
                    page_content.append(f"• {element_text}")
                elif element_type == 'Footer':
                    # Include footers as they may contain important info
                    page_content.append(f"\n*Footer: {element_text}*\n")
                elif element_type == 'Header':
                    # Include headers
                    page_content.append(f"\n*Header: {element_text}*\n")
                else:
                    # Include all other content types
                    page_content.append(element_text)
            
            if page_content:
                text_content.append(page_text + '\n'.join(page_content))
        
        # Combine all pages
        full_text = '\n\n'.join(text_content)
        
        # Add document metadata
        doc_header = f"# COMPLETE DOCUMENT EXTRACTION\n\n*Processed {len(page_elements)} pages with {len(elements)} total elements*\n\n"
        
        return doc_header + full_text
    
    def _fallback_pdf_processing(self, pdf_input, is_bytes: bool = False) -> str:
        """Fallback PDF processing using PyMuPDF if unstructured fails"""
        try:
            import fitz  # PyMuPDF
            
            logger.warning("Falling back to PyMuPDF for PDF processing")
            
            if is_bytes:
                doc = fitz.open(stream=pdf_input, filetype="pdf")
            else:
                doc = fitz.open(pdf_input)
            
            text_parts = []
            for page_num in range(len(doc)):
                page = doc[page_num]
                text = page.get_text()
                if text.strip():
                    text_parts.append(f"## Page {page_num + 1}\n\n{text.strip()}")
            
            doc.close()
            full_text = '\n\n'.join(text_parts)
            logger.info(f"Fallback PDF processing: {len(full_text)} characters")
            return full_text
            
        except Exception as e:
            logger.error(f"Fallback PDF processing also failed: {e}")
            return "Error: Could not process PDF content"
    
    def _chunk_text_by_structure(self, text: str, source: str) -> List[DocumentChunk]:
        """
        Chunk text based on structure while maintaining semantic coherence and completeness
        """
        chunks = []
        chunk_size = 3072  # Larger chunks for comprehensive documents like NIST
        chunk_overlap = 300  # More overlap to ensure continuity
        
        logger.info(f"Starting comprehensive chunking for {source}: {len(text):,} chars")
        
        # Handle both page-based and section-based splitting
        if "=== PAGE" in text:
            # Process page-by-page for comprehensive extraction
            pages = text.split("=== PAGE")
            logger.info(f"Found {len(pages)} pages in document")
            chunks.extend(self._chunk_by_pages(pages, source, chunk_size, chunk_overlap))
        else:
            # Fall back to section-based chunking
            chunks.extend(self._chunk_by_sections(text, source, chunk_size, chunk_overlap))
        
        # Validate chunking completeness
        total_chunk_chars = sum(len(chunk.text) for chunk in chunks)
        coverage_ratio = total_chunk_chars / len(text) if len(text) > 0 else 0
        
        logger.info(f"Chunking summary for {source}:")
        logger.info(f"  - Original: {len(text):,} chars")
        logger.info(f"  - Chunks: {len(chunks)} chunks")
        logger.info(f"  - Coverage: {coverage_ratio:.1%}")
        logger.info(f"  - Avg chunk size: {total_chunk_chars // len(chunks) if chunks else 0} chars")
        
        if coverage_ratio < 0.85:
            logger.warning(f"Low chunking coverage ({coverage_ratio:.1%}) for {source} - possible content loss")
        
        return chunks
    
    def _chunk_by_pages(self, pages, source: str, chunk_size: int, chunk_overlap: int) -> List[DocumentChunk]:
        """Chunk content by pages to ensure complete coverage"""
        chunks = []
        chunk_count = 0
        current_chunk = ""
        
        for i, page in enumerate(pages):
            if i == 0:  # Skip first empty split
                page_content = page
            else:
                page_content = f"=== PAGE{page}"
            
            page_content = page_content.strip()
            if not page_content:
                continue
            
            # Extract page number for metadata
            page_num_match = page_content.split('\n')[0] if page_content else ""
            page_num = i if i > 0 else 1
            
            # Try to add entire page to current chunk
            test_chunk = current_chunk + "\n\n" + page_content if current_chunk else page_content
            
            if len(test_chunk) <= chunk_size:
                current_chunk = test_chunk
            else:
                # Save current chunk if it has content
                if current_chunk.strip():
                    chunk_id = f"{source}_chunk_{chunk_count}"
                    chunks.append(DocumentChunk(
                        text=current_chunk.strip(),
                        source=source,
                        chunk_id=chunk_id,
                        metadata={
                            'chunk_index': chunk_count, 
                            'char_count': len(current_chunk),
                            'content_type': 'page_based'
                        }
                    ))
                    chunk_count += 1
                
                # Handle overlap
                if chunk_overlap > 0 and current_chunk:
                    overlap_text = current_chunk[-chunk_overlap:]
                    current_chunk = overlap_text + "\n\n" + page_content
                else:
                    current_chunk = page_content
                
                # If page is still too large, split it
                if len(current_chunk) > chunk_size:
                    page_chunks = self._split_large_content(current_chunk, source, chunk_count, chunk_size)
                    chunks.extend(page_chunks)
                    chunk_count += len(page_chunks)
                    current_chunk = ""
        
        # Add final chunk
        if current_chunk.strip():
            chunk_id = f"{source}_chunk_{chunk_count}"
            chunks.append(DocumentChunk(
                text=current_chunk.strip(),
                source=source,
                chunk_id=chunk_id,
                metadata={
                    'chunk_index': chunk_count, 
                    'char_count': len(current_chunk),
                    'content_type': 'page_based'
                }
            ))
        
        return chunks
    
    def _chunk_by_sections(self, text: str, source: str, chunk_size: int, chunk_overlap: int) -> List[DocumentChunk]:
        """Chunk content by sections (fallback method)"""
        chunks = []
        chunk_count = 0
        
        # Split by sections (marked by ##)
        sections = text.split('## ')
        current_chunk = ""
        
        for i, section in enumerate(sections):
            if i == 0 and not section.startswith('#'):
                section_text = section
            else:
                section_text = f"## {section}"
            
            section_text = section_text.strip()
            if not section_text:
                continue
            
            test_chunk = current_chunk + "\n\n" + section_text if current_chunk else section_text
            
            if len(test_chunk) <= chunk_size:
                current_chunk = test_chunk
            else:
                # Save current chunk
                if current_chunk.strip():
                    chunk_id = f"{source}_chunk_{chunk_count}"
                    chunks.append(DocumentChunk(
                        text=current_chunk.strip(),
                        source=source,
                        chunk_id=chunk_id,
                        metadata={
                            'chunk_index': chunk_count, 
                            'char_count': len(current_chunk),
                            'content_type': 'section_based'
                        }
                    ))
                    chunk_count += 1
                
                # Handle overlap
                if chunk_overlap > 0 and current_chunk:
                    overlap_text = current_chunk[-chunk_overlap:]
                    last_sentence = overlap_text.rfind('. ')
                    if last_sentence > 0:
                        overlap_text = overlap_text[last_sentence + 2:]
                    current_chunk = overlap_text + "\n\n" + section_text
                else:
                    current_chunk = section_text
                
                # Split large sections
                if len(current_chunk) > chunk_size:
                    section_chunks = self._split_large_content(current_chunk, source, chunk_count, chunk_size)
                    chunks.extend(section_chunks)
                    chunk_count += len(section_chunks)
                    current_chunk = ""
        
        # Add final chunk
        if current_chunk.strip():
            chunk_id = f"{source}_chunk_{chunk_count}"
            chunks.append(DocumentChunk(
                text=current_chunk.strip(),
                source=source,
                chunk_id=chunk_id,
                metadata={
                    'chunk_index': chunk_count, 
                    'char_count': len(current_chunk),
                    'content_type': 'section_based'
                }
            ))
        
        return chunks
    
    def _split_large_content(self, content: str, source: str, start_chunk_count: int, chunk_size: int) -> List[DocumentChunk]:
        """Split large content into smaller chunks while preserving context"""
        chunks = []
        paragraphs = content.split('\n\n')
        current_chunk = ""
        chunk_count = start_chunk_count
        
        for paragraph in paragraphs:
            test_chunk = current_chunk + "\n\n" + paragraph if current_chunk else paragraph
            
            if len(test_chunk) <= chunk_size:
                current_chunk = test_chunk
            else:
                # Save current chunk
                if current_chunk.strip():
                    chunk_id = f"{source}_chunk_{chunk_count}"
                    chunks.append(DocumentChunk(
                        text=current_chunk.strip(),
                        source=source,
                        chunk_id=chunk_id,
                        metadata={
                            'chunk_index': chunk_count, 
                            'char_count': len(current_chunk),
                            'content_type': 'split_large',
                            'parent_chunk': start_chunk_count
                        }
                    ))
                    chunk_count += 1
                
                current_chunk = paragraph
        
        # Add final chunk
        if current_chunk.strip():
            chunk_id = f"{source}_chunk_{chunk_count}"
            chunks.append(DocumentChunk(
                text=current_chunk.strip(),
                source=source,
                chunk_id=chunk_id,
                metadata={
                    'chunk_index': chunk_count, 
                    'char_count': len(current_chunk),
                    'content_type': 'split_large',
                    'parent_chunk': start_chunk_count
                }
            ))
        
        return chunks
    
    def add_document(self, text: str, source: str, metadata: Dict = None) -> int:
        """
        Add a document to the vector store
        
        Args:
            text: Document text content
            source: Document source/filename
            metadata: Optional metadata dictionary
            
        Returns:
            Number of chunks added
        """
        try:
            # Validate input
            if not text or not text.strip():
                logger.error(f"Document {source} is empty or contains no text")
                return 0
            
            # Generate document hash for deduplication
            doc_hash = hashlib.md5(f"{source}_{text}".encode()).hexdigest()
            
            # Check if document already exists
            if doc_hash in self.doc_metadata:
                logger.info(f"Document {source} already exists in vector store")
                return 0
            
            logger.info(f"Processing document: {source}")
            logger.info(f"   - Size: {len(text):,} characters ({len(text.split()):,} words)")
            
            # Chunk the document using structure-aware chunking
            chunks = self._chunk_text_by_structure(text, source)
            
            if not chunks:
                logger.error(f"No chunks generated for document: {source}")
                return 0
            
            logger.info(f"Generated {len(chunks)} chunks for {source}")
            
            # Generate embeddings for all chunks
            chunk_texts = [chunk.text for chunk in chunks]
            
            logger.info(f"Generating embeddings for {len(chunk_texts)} chunks...")
            try:
                embeddings = self.embedding_model.encode(chunk_texts, convert_to_numpy=True)
            except Exception as embed_error:
                logger.error(f"Failed to generate embeddings for {source}: {embed_error}")
                return 0
            
            # Add chunks to Qdrant
            points = []
            for i, (chunk, embedding) in enumerate(zip(chunks, embeddings)):
                point_id = self.next_point_id
                self.next_point_id += 1
                
                points.append(PointStruct(
                    id=point_id,
                    vector=embedding.tolist(),
                    payload={
                        "text": chunk.text,
                        "source": chunk.source,
                        "chunk_id": chunk.chunk_id,
                        "metadata": chunk.metadata,
                        "timestamp": chunk.timestamp
                    }
                ))
            
            # Insert points into Qdrant
            try:
                self.client.upsert(
                    collection_name=self.collection_name,
                    points=points
                )
                logger.info(f"Added {len(points)} points to Qdrant collection")
            except Exception as qdrant_error:
                logger.error(f"Failed to add points to Qdrant: {qdrant_error}")
                return 0
            
            # Store document chunks and metadata
            self.documents.extend(chunks)
            self.doc_metadata[doc_hash] = {
                'source': source,
                'chunk_count': len(chunks),
                'metadata': metadata or {},
                'timestamp': datetime.now().isoformat(),
                'original_size': len(text)
            }
            
            logger.info(f"Successfully added document {source} with {len(chunks)} chunks to vector store")
            return len(chunks)
            
        except Exception as e:
            logger.error(f"Critical error adding document {source} to vector store: {e}")
            import traceback
            logger.error(f"Full traceback: {traceback.format_exc()}")
            return 0
    
    def add_pdf_document(self, pdf_input, source: str, is_bytes: bool = False, metadata: Dict = None) -> int:
        """
        Add a PDF document to the vector store using unstructured.io
        
        Args:
            pdf_input: Path to PDF file or PDF bytes
            source: Document source/filename
            is_bytes: True if pdf_input is bytes, False if it's a file path
            metadata: Optional metadata dictionary
            
        Returns:
            Number of chunks added
        """
        try:
            logger.info(f"Processing PDF document: {source} with unstructured.io")
            
            # Process PDF with unstructured.io
            text_content = self._process_pdf_with_unstructured(pdf_input, is_bytes)
            
            if not text_content or text_content.startswith("Error:"):
                logger.error(f"Failed to extract content from PDF: {source}")
                return 0
            
            # Add the processed text as a regular document
            pdf_metadata = metadata or {}
            pdf_metadata['document_type'] = 'pdf'
            pdf_metadata['processing_method'] = 'unstructured.io'
            
            return self.add_document(text_content, source, pdf_metadata)
            
        except Exception as e:
            logger.error(f"Error adding PDF document {source}: {e}")
            import traceback
            logger.error(f"Full traceback: {traceback.format_exc()}")
            return 0
    
    def search(self, query: str, top_k: int = 5, score_threshold: float = 0.3) -> List[Tuple[DocumentChunk, float]]:
        """
        Search for relevant document chunks
        
        Args:
            query: Search query
            top_k: Number of top results to return
            score_threshold: Minimum similarity score
            
        Returns:
            List of (DocumentChunk, score) tuples
        """
        try:
            # Check if collection has any points
            collection_info = self.client.get_collection(self.collection_name)
            if collection_info.points_count == 0:
                logger.info("Vector store is empty, no search results")
                return []
            
            # Generate query embedding
            query_embedding = self.embedding_model.encode([query], convert_to_numpy=True)[0]
            
            # Search Qdrant
            search_results = self.client.search(
                collection_name=self.collection_name,
                query_vector=query_embedding.tolist(),
                limit=top_k,
                score_threshold=score_threshold
            )
            
            # Convert results to DocumentChunk objects
            results = []
            for result in search_results:
                payload = result.payload
                chunk = DocumentChunk(
                    text=payload["text"],
                    source=payload["source"],
                    chunk_id=payload["chunk_id"],
                    metadata=payload.get("metadata", {})
                )
                results.append((chunk, result.score))
            
            logger.info(f"Vector search for '{query[:50]}...' returned {len(results)} results")
            return results
            
        except Exception as e:
            logger.error(f"Error searching vector store: {e}")
            return []
    
    def get_all_documents_context(self, max_chars: int = 50000) -> str:
        """
        Get context from all documents in the store
        
        Args:
            max_chars: Maximum total characters in context
            
        Returns:
            Formatted context string with excerpts from all documents
        """
        try:
            if not self.documents:
                return ""
            
            context_parts = []
            total_chars = 0
            
            # Group documents by source
            by_source = {}
            for doc in self.documents:
                if doc.source not in by_source:
                    by_source[doc.source] = []
                by_source[doc.source].append(doc)
            
            logger.info(f"Processing {len(by_source)} documents for comprehensive analysis")
            
            # Process all chunks from each document
            for source, chunks in by_source.items():
                source_text = f"\n--- From {source} ---\n"
                context_parts.append(source_text)
                total_chars += len(source_text)
                
                for i, chunk in enumerate(chunks):
                    chunk_text = chunk.text + "\n\n"
                    
                    if total_chars + len(chunk_text) > max_chars:
                        remaining_chars = max_chars - total_chars
                        if remaining_chars > 100:
                            chunk_text = chunk.text[:remaining_chars - 50] + "...\n\n"
                            context_parts.append(chunk_text)
                            total_chars += len(chunk_text)
                        break
                    
                    context_parts.append(chunk_text)
                    total_chars += len(chunk_text)
                
                if total_chars >= max_chars:
                    break
            
            context = "".join(context_parts).strip()
            logger.info(f"Generated comprehensive context: {len(context)} chars from {len(by_source)} documents")
            return context
            
        except Exception as e:
            logger.error(f"Error generating comprehensive context: {e}")
            return ""
    
    def get_context_for_query(self, query: str, max_chunks: int = 50, max_chars: int = 50000) -> str:
        """
        Get relevant context for a query, formatted for LLM consumption
        
        Args:
            query: User query
            max_chunks: Maximum number of chunks to include
            max_chars: Maximum total characters in context
            
        Returns:
            Formatted context string
        """
        try:
            results = self.search(query, top_k=max_chunks)
            
            if not results:
                return ""
            
            context_parts = []
            total_chars = 0
            
            # Group results by source document
            by_source = {}
            for chunk, score in results:
                source = chunk.source
                if source not in by_source:
                    by_source[source] = []
                by_source[source].append((chunk, score))
            
            # Format context by source
            for source, source_chunks in by_source.items():
                source_text = f"\n--- From {source} ---\n"
                
                for chunk, score in source_chunks:
                    chunk_text = f"{chunk.text}\n"
                    
                    # Check if adding this chunk would exceed character limit
                    if total_chars + len(source_text) + len(chunk_text) > max_chars:
                        break
                    
                    if not any(source in part for part in context_parts):
                        context_parts.append(source_text)
                        total_chars += len(source_text)
                    
                    context_parts.append(chunk_text)
                    total_chars += len(chunk_text)
                
                if total_chars >= max_chars:
                    break
            
            context = "".join(context_parts).strip()
            logger.info(f"Generated context: {len(context)} chars from {len(by_source)} documents")
            return context
            
        except Exception as e:
            logger.error(f"Error generating context for query: {e}")
            return ""
    
    def get_full_documents_context(self, max_chars: int = 200000) -> str:
        """
        Get ALL content from ALL documents with minimal limits for maximum comprehensiveness
        
        Args:
            max_chars: Very high limit for comprehensive analysis
            
        Returns:
            Complete document content formatted for LLM consumption
        """
        try:
            if not self.documents:
                return ""
            
            context_parts = []
            total_chars = 0
            
            # Group documents by source
            by_source = {}
            for doc in self.documents:
                if doc.source not in by_source:
                    by_source[doc.source] = []
                by_source[doc.source].append(doc)
            
            logger.info(f"Retrieving COMPLETE content from {len(by_source)} documents (limit: {max_chars:,} chars)")
            
            # Process ALL chunks from ALL documents
            for source, chunks in by_source.items():
                source_text = f"\n=== COMPLETE CONTENT FROM {source.upper()} ===\n\n"
                context_parts.append(source_text)
                total_chars += len(source_text)
                
                logger.info(f"Adding ALL {len(chunks)} chunks from {source}")
                
                chunks_added = 0
                for chunk in chunks:
                    chunk_text = chunk.text + "\n\n"
                    
                    # Memory safety check
                    if len(chunk_text) > 50000:
                        logger.warning(f"Skipping very large chunk ({len(chunk_text):,} chars) from {source}")
                        continue
                    
                    if total_chars + len(chunk_text) > max_chars:
                        logger.warning(f"Reached {max_chars:,} char limit after {chunks_added}/{len(chunks)} chunks from {source}")
                        break
                    
                    context_parts.append(chunk_text)
                    total_chars += len(chunk_text)
                    chunks_added += 1
                
                logger.info(f"Added {chunks_added}/{len(chunks)} chunks from {source}")
                
                # Add document separator
                if total_chars < max_chars - 100:
                    separator = f"\n--- End of {source} ---\n\n"
                    context_parts.append(separator)
                    total_chars += len(separator)
            
            context = "".join(context_parts).strip()
            logger.info(f"Generated COMPLETE context: {len(context):,} chars from {len(by_source)} documents")
            return context
            
        except Exception as e:
            logger.error(f"Error generating complete documents context: {e}")
            return ""
    
    def remove_document(self, source: str) -> bool:
        """
        Remove a document from the vector store
        """
        try:
            # Find and remove points with matching source
            search_result = self.client.scroll(
                collection_name=self.collection_name,
                scroll_filter=models.Filter(
                    must=[
                        models.FieldCondition(
                            key="source",
                            match=models.MatchValue(value=source)
                        )
                    ]
                ),
                limit=10000  # Large limit to get all matching points
            )
            
            point_ids = [point.id for point in search_result[0]]
            
            if point_ids:
                self.client.delete(
                    collection_name=self.collection_name,
                    points_selector=models.PointIdsList(points=point_ids)
                )
                logger.info(f"Removed {len(point_ids)} points for document {source}")
            
            # Remove from local tracking
            self.documents = [doc for doc in self.documents if doc.source != source]
            
            # Remove from metadata
            doc_hash_to_remove = None
            for doc_hash, metadata in self.doc_metadata.items():
                if metadata['source'] == source:
                    doc_hash_to_remove = doc_hash
                    break
            
            if doc_hash_to_remove:
                del self.doc_metadata[doc_hash_to_remove]
            
            logger.info(f"Removed document {source} from vector store")
            return True
            
        except Exception as e:
            logger.error(f"Error removing document {source}: {e}")
            return False
    
    def get_stats(self) -> Dict[str, Any]:
        """Get vector store statistics"""
        try:
            collection_info = self.client.get_collection(self.collection_name)
            
            active_docs = []
            total_chunks = 0
            
            for doc_hash, meta in self.doc_metadata.items():
                doc_info = {
                    'source': meta['source'],
                    'chunks': meta['chunk_count'],
                    'original_size': meta.get('original_size', 0),
                    'timestamp': meta['timestamp']
                }
                active_docs.append(doc_info)
                total_chunks += meta['chunk_count']
            
            return {
                'session_id': self.session_id,
                'total_documents': len(active_docs),
                'total_chunks': total_chunks,
                'index_size': collection_info.points_count,
                'embedding_model': self.embedding_model_name,
                'embedding_dim': self.embedding_dim,
                'documents': active_docs
            }
        except Exception as e:
            logger.error(f"Error getting stats: {e}")
            return {'error': str(e)}
    
    def get_document_sample(self, source: str, max_chunks: int = 5) -> List[str]:
        """Get sample chunks from a specific document for debugging"""
        try:
            sample_chunks = []
            for chunk in self.documents:
                if chunk.source == source:
                    sample_chunks.append({
                        'chunk_id': chunk.chunk_id,
                        'text_preview': chunk.text[:200] + "..." if len(chunk.text) > 200 else chunk.text,
                        'text_length': len(chunk.text)
                    })
                    if len(sample_chunks) >= max_chunks:
                        break
            
            logger.info(f"Found {len(sample_chunks)} matching chunks for '{source}'")
            return sample_chunks
        except Exception as e:
            logger.error(f"Error getting document sample for {source}: {e}")
            return []
    
    def clear(self):
        """Clear all documents from the vector store"""
        try:
            # Delete the collection and recreate it
            self.client.delete_collection(self.collection_name)
            self.client.create_collection(
                collection_name=self.collection_name,
                vectors_config=VectorParams(
                    size=self.embedding_dim,
                    distance=Distance.COSINE
                )
            )
            
            # Clear local tracking
            self.documents = []
            self.doc_metadata = {}
            self.next_point_id = 0
            
            logger.info("Cleared vector store")
        except Exception as e:
            logger.error(f"Error clearing vector store: {e}")

class VectorStoreManager:
    """Manages vector store instances per session"""
    
    def __init__(self):
        """Initialize vector store manager"""
        self.sessions: Dict[str, QdrantVectorStore] = {}
        logger.info("Qdrant Vector store manager initialized")
    
    def get_or_create_store(self, session_id: str) -> Optional[QdrantVectorStore]:
        """Get or create a vector store for a session"""
        try:
            if not VECTOR_STORE_AVAILABLE:
                logger.warning("Qdrant not available, cannot create vector store")
                return None
            
            if session_id not in self.sessions:
                self.sessions[session_id] = QdrantVectorStore(session_id)
                logger.info(f"Created new Qdrant vector store for session: {session_id}")
            
            return self.sessions[session_id]
            
        except Exception as e:
            logger.error(f"Error creating vector store for session {session_id}: {e}")
            return None
    
    def add_document_from_content(self, session_id: str, content: str, source_filename: str, metadata: dict = None) -> bool:
        """Add a document to vector store from content string"""
        try:
            store = self.get_or_create_store(session_id)
            if not store:
                logger.error(f"Failed to get/create vector store for session {session_id}")
                return False
            
            # Check if this is PDF content (markdown converted from PDF)
            if metadata and metadata.get('type') == 'pdf_markdown':
                # This is PDF content converted to markdown, add as regular document
                chunks_added = store.add_document(content, source_filename, metadata or {})
            else:
                # Regular text content
                chunks_added = store.add_document(content, source_filename, metadata or {})
            
            logger.info(f"Added document {source_filename} to vector store: {chunks_added} chunks")
            return chunks_added > 0
            
        except Exception as e:
            logger.error(f"Error adding document {source_filename} to vector store: {e}")
            import traceback
            logger.error(f"Full traceback: {traceback.format_exc()}")
            return False
    
    def add_pdf_document(self, session_id: str, pdf_input, source_filename: str, is_bytes: bool = False, metadata: dict = None) -> bool:
        """Add a PDF document to vector store using unstructured.io processing"""
        try:
            store = self.get_or_create_store(session_id)
            if not store:
                logger.error(f"Failed to get/create vector store for session {session_id}")
                return False
            
            chunks_added = store.add_pdf_document(pdf_input, source_filename, is_bytes, metadata or {})
            logger.info(f"Added PDF document {source_filename} to vector store: {chunks_added} chunks")
            return chunks_added > 0
            
        except Exception as e:
            logger.error(f"Error adding PDF document {source_filename} to vector store: {e}")
            import traceback
            logger.error(f"Full traceback: {traceback.format_exc()}")
            return False
    
    def add_pdf_document_from_s3(self, session_id: str, s3_key: str, s3_bucket: str, source_filename: str, metadata: dict = None) -> bool:
        """Add a PDF document to vector store from S3 using unstructured.io processing"""
        try:
            store = self.get_or_create_store(session_id)
            if not store:
                logger.error(f"Failed to get/create vector store for session {session_id}")
                return False
            
            # Download PDF from S3 first
            try:
                import boto3
                s3_client = boto3.client('s3')
                response = s3_client.get_object(Bucket=s3_bucket, Key=s3_key)
                pdf_bytes = response['Body'].read()
                logger.info(f"Downloaded PDF from S3: {source_filename} ({len(pdf_bytes)} bytes)")
                
                # Process the PDF bytes with unstructured.io
                chunks_added = store.add_pdf_document(pdf_bytes, source_filename, is_bytes=True, metadata=metadata or {})
                logger.info(f"Added PDF document {source_filename} from S3 to vector store: {chunks_added} chunks")
                return chunks_added > 0
                
            except Exception as s3_error:
                logger.error(f"Error downloading PDF from S3 {s3_bucket}/{s3_key}: {s3_error}")
                return False
            
        except Exception as e:
            logger.error(f"Error adding PDF document {source_filename} from S3 to vector store: {e}")
            import traceback
            logger.error(f"Full traceback: {traceback.format_exc()}")
            return False
    
    def add_document_from_file(self, session_id: str, file_path: str, source_filename: str) -> bool:
        """Add a document to vector store from local file"""
        try:
            store = self.get_or_create_store(session_id)
            if not store:
                return False
            
            # Check if it's a PDF file
            if file_path.lower().endswith('.pdf'):
                chunks_added = store.add_pdf_document(file_path, source_filename, is_bytes=False, metadata={'file_path': file_path})
            else:
                # Read text file
                with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                    content = f.read()
                chunks_added = store.add_document(content, source_filename, {'file_path': file_path})
            
            logger.info(f"Added document {source_filename} from file to vector store: {chunks_added} chunks")
            return chunks_added > 0
            
        except Exception as e:
            logger.error(f"Error adding document from file to vector store: {e}")
            return False
    
    def search_session(self, session_id: str, query: str, top_k: int = 5) -> List[Tuple[DocumentChunk, float]]:
        """Search documents in a session's vector store"""
        try:
            if session_id not in self.sessions:
                logger.info(f"No vector store found for session: {session_id}")
                return []
            
            return self.sessions[session_id].search(query, top_k)
            
        except Exception as e:
            logger.error(f"Error searching session {session_id}: {e}")
            return []
    
    def get_all_documents_for_session(self, session_id: str, max_chars: int = 50000, force_full_retrieval: bool = False) -> str:
        """Get context from all documents in a session's vector store"""
        try:
            if session_id not in self.sessions:
                logger.info(f"No vector store found for session: {session_id}")
                return ""
            
            if force_full_retrieval:
                return self.sessions[session_id].get_full_documents_context(max_chars)
            else:
                return self.sessions[session_id].get_all_documents_context(max_chars)
            
        except Exception as e:
            logger.error(f"Error getting all documents for session {session_id}: {e}")
            return ""

    def get_context_for_session(self, session_id: str, query: str, max_chunks: int = 20, max_chars: int = 25000) -> str:
        """Get relevant context for a query from session's vector store"""
        try:
            if session_id not in self.sessions:
                logger.info(f"No vector store found for session: {session_id}")
                return ""
            
            return self.sessions[session_id].get_context_for_query(query, max_chunks, max_chars)
            
        except Exception as e:
            logger.error(f"Error getting context for session {session_id}: {e}")
            return ""
    
    def remove_document_from_session(self, session_id: str, source_filename: str) -> bool:
        """Remove a document from session's vector store"""
        try:
            if session_id not in self.sessions:
                logger.info(f"No vector store found for session: {session_id}")
                return False
            
            return self.sessions[session_id].remove_document(source_filename)
            
        except Exception as e:
            logger.error(f"Error removing document from session {session_id}: {e}")
            return False
    
    def clear_session(self, session_id: str):
        """Clear a session's vector store"""
        try:
            if session_id in self.sessions:
                self.sessions[session_id].clear()
                del self.sessions[session_id]
                logger.info(f"Cleared vector store for session: {session_id}")
        except Exception as e:
            logger.error(f"Error clearing session {session_id}: {e}")
    
    def get_session_stats(self, session_id: str) -> Dict[str, Any]:
        """Get statistics for a session's vector store"""
        try:
            if session_id not in self.sessions:
                return {'error': 'Session not found'}
            
            return self.sessions[session_id].get_stats()
            
        except Exception as e:
            logger.error(f"Error getting stats for session {session_id}: {e}")
            return {'error': str(e)}
    
    def get_document_sample(self, session_id: str, source: str, max_chunks: int = 5) -> List[str]:
        """Get sample chunks from a specific document for debugging"""
        try:
            if session_id not in self.sessions:
                return []
            
            return self.sessions[session_id].get_document_sample(source, max_chunks)
            
        except Exception as e:
            logger.error(f"Error getting document sample for session {session_id}: {e}")
            return []

# Global vector store manager instance
vector_store_manager = None

def initialize_vector_store_manager() -> Optional[VectorStoreManager]:
    """Initialize the global vector store manager"""
    global vector_store_manager
    
    try:
        if not VECTOR_STORE_AVAILABLE:
            logger.warning("Qdrant/Unstructured not available, vector store disabled")
            return None
        
        vector_store_manager = VectorStoreManager()
        logger.info("Initialized Qdrant vector store manager")
        return vector_store_manager
        
    except Exception as e:
        logger.error(f"Error initializing vector store manager: {e}")
        return None

def get_vector_store_manager() -> Optional[VectorStoreManager]:
    """Get the global vector store manager instance"""
    return vector_store_manager

# Utility function to check if vector store is available
def is_vector_store_available() -> bool:
    """Check if vector store functionality is available"""
    return VECTOR_STORE_AVAILABLE and vector_store_manager is not None