"""
FAISS Vector Store for BedBot - Document Embeddings and Retrieval
"""

import os
import json
import logging
import tempfile
import pickle
import hashlib
from datetime import datetime
from typing import List, Dict, Tuple, Optional, Any

# Set up logger
logger = logging.getLogger(__name__)

try:
    import faiss
    import numpy as np
    from sentence_transformers import SentenceTransformer
    FAISS_AVAILABLE = True
    logger.info("FAISS and SentenceTransformers loaded successfully")
except ImportError as e:
    FAISS_AVAILABLE = False
    logger.warning(f"FAISS dependencies not available: {e}")
    logger.warning("To enable vector store: pip install faiss-cpu sentence-transformers")

class DocumentChunk:
    """Represents a chunk of text from a document"""
    def __init__(self, text: str, source: str, chunk_id: str, metadata: Dict = None):
        self.text = text
        self.source = source  # filename
        self.chunk_id = chunk_id
        self.metadata = metadata or {}
        self.timestamp = datetime.now().isoformat()

class FAISSVectorStore:
    """FAISS-based vector store for document embeddings and retrieval"""
    
    def __init__(self, session_id: str, embedding_model: str = "all-MiniLM-L6-v2", 
                 chunk_size: int = 2048, chunk_overlap: int = 200):
        """
        Initialize FAISS vector store
        
        Args:
            session_id: Unique session identifier
            embedding_model: SentenceTransformer model name
            chunk_size: Maximum characters per chunk
            chunk_overlap: Character overlap between chunks
        """
        if not FAISS_AVAILABLE:
            raise ImportError("FAISS dependencies not available. Install with: pip install faiss-cpu sentence-transformers")
        
        self.session_id = session_id
        self.embedding_model_name = embedding_model
        self.chunk_size = chunk_size
        self.chunk_overlap = chunk_overlap
        
        # Initialize embedding model
        try:
            self.embedding_model = SentenceTransformer(embedding_model)
            self.embedding_dim = self.embedding_model.get_sentence_embedding_dimension()
            logger.info(f"Loaded embedding model: {embedding_model} (dim: {self.embedding_dim})")
        except Exception as e:
            logger.error(f"Failed to load embedding model {embedding_model}: {e}")
            raise
        
        # Initialize FAISS index
        self.index = faiss.IndexFlatIP(self.embedding_dim)  # Inner product (cosine similarity)
        self.documents: List[DocumentChunk] = []
        self.doc_metadata: Dict[str, Any] = {}
        
        # State tracking
        self.is_dirty = False  # Track if index needs saving
        
        logger.info(f"Initialized FAISS vector store for session: {session_id}")
    
    def _chunk_text(self, text: str, source: str) -> List[DocumentChunk]:
        """Split text into overlapping chunks optimized for large documents and context windows"""
        chunks = []
        
        # For large documents, use paragraph-aware chunking for better semantic coherence
        paragraphs = text.split('\n\n')
        current_chunk = ""
        chunk_count = 0
        
        for paragraph in paragraphs:
            # Clean paragraph
            paragraph = paragraph.strip()
            if not paragraph:
                continue
                
            # Try adding paragraph to current chunk
            test_chunk = current_chunk + "\n\n" + paragraph if current_chunk else paragraph
            
            if len(test_chunk) <= self.chunk_size:
                current_chunk = test_chunk
            else:
                # Current chunk is getting too large
                if current_chunk.strip():
                    # Save current chunk
                    chunk_id = f"{source}_{chunk_count}"
                    chunks.append(DocumentChunk(
                        text=current_chunk.strip(),
                        source=source,
                        chunk_id=chunk_id,
                        metadata={'chunk_index': chunk_count, 'char_count': len(current_chunk)}
                    ))
                    chunk_count += 1
                
                # Handle overlap for continuity
                if self.chunk_overlap > 0 and current_chunk:
                    # Take last overlap characters, but try to break at sentence boundary
                    overlap_text = current_chunk[-self.chunk_overlap:]
                    # Find last sentence boundary in overlap
                    last_sentence = overlap_text.rfind('. ')
                    if last_sentence > 0:
                        overlap_text = overlap_text[last_sentence + 2:]
                    current_chunk = overlap_text + "\n\n" + paragraph
                else:
                    current_chunk = paragraph
                
                # If single paragraph is larger than chunk_size, split by sentences
                if len(current_chunk) > self.chunk_size:
                    sentences = current_chunk.split('. ')
                    sentence_chunk = ""
                    
                    for sentence in sentences:
                        test_sentence_chunk = sentence_chunk + sentence + ". "
                        if len(test_sentence_chunk) <= self.chunk_size:
                            sentence_chunk = test_sentence_chunk
                        else:
                            # Save sentence chunk
                            if sentence_chunk.strip():
                                chunk_id = f"{source}_{chunk_count}"
                                chunks.append(DocumentChunk(
                                    text=sentence_chunk.strip(),
                                    source=source,
                                    chunk_id=chunk_id,
                                    metadata={'chunk_index': chunk_count, 'char_count': len(sentence_chunk)}
                                ))
                                chunk_count += 1
                            sentence_chunk = sentence + ". "
                    
                    current_chunk = sentence_chunk
        
        # Add final chunk if not empty
        if current_chunk.strip():
            chunk_id = f"{source}_{chunk_count}"
            chunks.append(DocumentChunk(
                text=current_chunk.strip(),
                source=source,
                chunk_id=chunk_id,
                metadata={'chunk_index': chunk_count, 'char_count': len(current_chunk)}
            ))
        
        logger.info(f"Chunked document {source} ({len(text)} chars) into {len(chunks)} chunks")
        
        # Debug: log first chunk preview for NIST documents
        if chunks and ("nist" in source.lower() or "800-171" in source.lower()):
            first_chunk_preview = chunks[0].text[:200] + "..." if len(chunks[0].text) > 200 else chunks[0].text
            logger.info(f"🔍 First chunk preview for {source}: {first_chunk_preview}")
            
            # Check if we're getting substantive content or just headers
            if "cover" in first_chunk_preview.lower() or "publication" in first_chunk_preview.lower():
                logger.warning(f"⚠️ {source} first chunk appears to be header/cover content - check PDF conversion")
        
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
                logger.error(f"❌ Document {source} is empty or contains no text")
                return 0
            
            if len(text) < 50:  # Very short documents might indicate processing issues
                logger.warning(f"⚠️ Document {source} is very short ({len(text)} chars) - possible processing issue")
            
            # Generate document hash for deduplication
            doc_hash = hashlib.md5(f"{source}_{text}".encode()).hexdigest()
            
            # Check if document already exists
            if doc_hash in self.doc_metadata:
                logger.info(f"Document {source} already exists in vector store")
                return 0
            
            # Log document stats
            logger.info(f"📄 Processing document: {source}")
            logger.info(f"   - Size: {len(text):,} characters ({len(text.split()):,} words)")
            logger.info(f"   - Expected chunks: ~{len(text) // self.chunk_size}")
            
            # Chunk the document
            chunks = self._chunk_text(text, source)
            
            if not chunks:
                logger.error(f"❌ No chunks generated for document: {source}")
                return 0
            
            # Validate chunk generation
            total_chunk_chars = sum(len(chunk.text) for chunk in chunks)
            coverage_ratio = total_chunk_chars / len(text)
            
            if coverage_ratio < 0.8:  # Less than 80% coverage suggests issues
                logger.warning(f"⚠️ Document {source} chunking coverage is low: {coverage_ratio:.1%}")
                logger.warning(f"   - Original: {len(text):,} chars")
                logger.warning(f"   - Chunks total: {total_chunk_chars:,} chars")
            
            logger.info(f"✅ Chunked document {source}: {len(chunks)} chunks, {coverage_ratio:.1%} coverage")
            
            # Generate embeddings for all chunks with memory safety
            chunk_texts = [chunk.text for chunk in chunks]
            
            logger.info(f"🔄 Generating embeddings for {len(chunk_texts)} chunks...")
            try:
                # Process in smaller batches to prevent memory issues
                batch_size = 50  # Process 50 chunks at a time
                all_embeddings = []
                
                for i in range(0, len(chunk_texts), batch_size):
                    batch = chunk_texts[i:i + batch_size]
                    logger.info(f"Processing embedding batch {i//batch_size + 1}/{(len(chunk_texts) + batch_size - 1)//batch_size}")
                    
                    batch_embeddings = self.embedding_model.encode(batch, convert_to_numpy=True)
                    all_embeddings.append(batch_embeddings)
                
                # Combine all embeddings
                import numpy as np
                embeddings = np.vstack(all_embeddings) if all_embeddings else np.array([])
                
            except Exception as embed_error:
                logger.error(f"❌ Failed to generate embeddings for {source}: {embed_error}")
                return 0
            
            # Validate embeddings
            if embeddings.shape[0] != len(chunks):
                logger.error(f"❌ Embedding count mismatch for {source}: {embeddings.shape[0]} != {len(chunks)}")
                return 0
            
            # Normalize embeddings for cosine similarity
            faiss.normalize_L2(embeddings)
            
            # Add to FAISS index
            try:
                self.index.add(embeddings)
            except Exception as faiss_error:
                logger.error(f"❌ Failed to add embeddings to FAISS index for {source}: {faiss_error}")
                return 0
            
            # Store document chunks
            start_idx = len(self.documents)
            self.documents.extend(chunks)
            
            # Store document metadata
            self.doc_metadata[doc_hash] = {
                'source': source,
                'chunk_count': len(chunks),
                'start_idx': start_idx,
                'end_idx': len(self.documents),
                'metadata': metadata or {},
                'timestamp': datetime.now().isoformat(),
                'original_size': len(text),
                'coverage_ratio': coverage_ratio
            }
            
            self.is_dirty = True
            logger.info(f"✅ Successfully added document {source} with {len(chunks)} chunks to vector store")
            return len(chunks)
            
        except Exception as e:
            logger.error(f"❌ Critical error adding document {source} to vector store: {e}")
            import traceback
            logger.error(f"Full traceback: {traceback.format_exc()}")
            return 0
    
    def search(self, query: str, top_k: int = 5, score_threshold: float = 0.1) -> List[Tuple[DocumentChunk, float]]:
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
            if self.index.ntotal == 0:
                logger.info("Vector store is empty, no search results")
                return []
            
            # Generate query embedding
            query_embedding = self.embedding_model.encode([query], convert_to_numpy=True)
            faiss.normalize_L2(query_embedding)
            
            # Search FAISS index
            scores, indices = self.index.search(query_embedding, min(top_k, self.index.ntotal))
            
            # Filter and format results
            results = []
            for score, idx in zip(scores[0], indices[0]):
                if idx >= 0 and score >= score_threshold and idx < len(self.documents):
                    results.append((self.documents[idx], float(score)))
            
            logger.info(f"Vector search for '{query[:50]}...' returned {len(results)} results")
            return results
            
        except Exception as e:
            logger.error(f"Error searching vector store: {e}")
            return []
    
    def get_all_documents_context(self, max_chars: int = 50000) -> str:
        """
        Get context from all documents in the store (for comprehensive analysis)
        
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
            
            num_documents = len(by_source)
            logger.info(f"Processing {num_documents} documents for comprehensive analysis - using ALL available chunks")
            logger.info(f"Maximum context limit: {max_chars:,} characters")
            
            # Use ALL chunks from each document for comprehensive coverage
            for source, chunks in by_source.items():
                source_text = f"\n--- From {source} ---\n"
                context_parts.append(source_text)
                total_chars += len(source_text)
                
                logger.info(f"Processing ALL {len(chunks)} chunks from {source}")
                
                # Use ALL chunks from each document
                for i, chunk in enumerate(chunks):
                    # Add full chunk text without arbitrary limits
                    chunk_text = chunk.text + "\n\n"
                    
                    # Check if adding this chunk would exceed the total limit
                    if total_chars + len(chunk_text) > max_chars:
                        # Try to fit as much as possible
                        remaining_chars = max_chars - total_chars
                        if remaining_chars > 100:  # Only add if we have meaningful space
                            chunk_text = chunk.text[:remaining_chars - 50] + "...\n\n"
                            context_parts.append(chunk_text)
                            total_chars += len(chunk_text)
                        logger.info(f"Reached character limit after {i+1}/{len(chunks)} chunks from {source}")
                        break
                    
                    context_parts.append(chunk_text)
                    total_chars += len(chunk_text)
                
                # Log completion for this document
                logger.info(f"Added complete content from {source}: {len([c for c in chunks if any(c.text in part for part in context_parts)])} chunks included")
                
                # Stop if we've reached the total limit
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
            max_chars: Very high limit for comprehensive analysis (default 500K)
            
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
            
            logger.info(f"🔄 Retrieving COMPLETE content from {len(by_source)} documents (limit: {max_chars:,} chars)")
            
            # Process ALL chunks from ALL documents
            for source, chunks in by_source.items():
                source_text = f"\n=== COMPLETE CONTENT FROM {source.upper()} ===\n\n"
                context_parts.append(source_text)
                total_chars += len(source_text)
                
                logger.info(f"📋 Adding ALL {len(chunks)} chunks from {source}")
                
                chunks_added = 0
                for chunk in chunks:
                    chunk_text = chunk.text + "\n\n"
                    
                    # Memory safety: check both total limit and individual chunk size
                    if len(chunk_text) > 50000:  # Skip extremely large chunks that might cause issues
                        logger.warning(f"⚠️ Skipping very large chunk ({len(chunk_text):,} chars) from {source}")
                        continue
                    
                    # Only stop if we absolutely must (very high limit)
                    if total_chars + len(chunk_text) > max_chars:
                        logger.warning(f"⚠️ Reached {max_chars:,} char limit after {chunks_added}/{len(chunks)} chunks from {source}")
                        break
                    
                    context_parts.append(chunk_text)
                    total_chars += len(chunk_text)
                    chunks_added += 1
                    
                    # Safety check: prevent memory issues with too many parts
                    if len(context_parts) > 1000:  # Reasonable limit on list size
                        logger.warning(f"⚠️ Reached maximum context parts limit (1000) for {source}")
                        break
                
                logger.info(f"✅ Added {chunks_added}/{len(chunks)} chunks from {source} ({sum(len(c.text) for c in chunks[:chunks_added]):,} chars)")
                
                # Add document separator
                if total_chars < max_chars - 100:
                    separator = f"\n--- End of {source} ---\n\n"
                    context_parts.append(separator)
                    total_chars += len(separator)
            
            context = "".join(context_parts).strip()
            logger.info(f"🎯 Generated COMPLETE context: {len(context):,} chars from {len(by_source)} documents")
            return context
            
        except Exception as e:
            logger.error(f"❌ Error generating complete documents context: {e}")
            return ""
    
    def remove_document(self, source: str) -> bool:
        """
        Remove a document from the vector store
        Note: FAISS doesn't support deletion, so this marks as removed
        """
        try:
            # Find document by source
            doc_hash = None
            for hash_key, metadata in self.doc_metadata.items():
                if metadata['source'] == source:
                    doc_hash = hash_key
                    break
            
            if not doc_hash:
                logger.info(f"Document {source} not found in vector store")
                return False
            
            # Mark as removed (FAISS doesn't support true deletion)
            self.doc_metadata[doc_hash]['removed'] = True
            self.is_dirty = True
            
            logger.info(f"Marked document {source} as removed from vector store")
            return True
            
        except Exception as e:
            logger.error(f"Error removing document {source}: {e}")
            return False
    
    def get_stats(self) -> Dict[str, Any]:
        """Get vector store statistics with detailed document info"""
        active_docs = []
        total_chunks = 0
        
        for doc_hash, meta in self.doc_metadata.items():
            if not meta.get('removed', False):
                doc_info = {
                    'source': meta['source'],
                    'chunks': meta['chunk_count'],
                    'original_size': meta.get('original_size', 0),
                    'coverage': meta.get('coverage_ratio', 0),
                    'timestamp': meta['timestamp']
                }
                active_docs.append(doc_info)
                total_chunks += meta['chunk_count']
        
        return {
            'session_id': self.session_id,
            'total_documents': len(active_docs),
            'total_chunks': total_chunks,
            'index_size': self.index.ntotal,
            'embedding_model': self.embedding_model_name,
            'embedding_dim': self.embedding_dim,
            'is_dirty': self.is_dirty,
            'documents': active_docs  # Detailed document list
        }
    
    def get_document_sample(self, source: str, max_chunks: int = 5) -> List[str]:
        """Get sample chunks from a specific document for debugging"""
        try:
            logger.info(f"🔍 Searching for chunks from source: '{source}'")
            logger.info(f"🔍 Total documents in store: {len(self.documents)}")
            
            # Debug: show all unique sources
            all_sources = set(chunk.source for chunk in self.documents)
            logger.info(f"🔍 All sources in vector store: {list(all_sources)}")
            
            sample_chunks = []
            for chunk in self.documents:
                logger.debug(f"🔍 Comparing chunk source '{chunk.source}' with target '{source}'")
                if chunk.source == source:
                    sample_chunks.append({
                        'chunk_id': chunk.chunk_id,
                        'text_preview': chunk.text[:200] + "..." if len(chunk.text) > 200 else chunk.text,
                        'text_length': len(chunk.text)
                    })
                    if len(sample_chunks) >= max_chunks:
                        break
            
            logger.info(f"🔍 Found {len(sample_chunks)} matching chunks for '{source}'")
            return sample_chunks
        except Exception as e:
            logger.error(f"❌ Error getting document sample for {source}: {e}")
            import traceback
            logger.error(f"Full traceback: {traceback.format_exc()}")
            return []
    
    def clear(self):
        """Clear all documents from the vector store"""
        try:
            self.index = faiss.IndexFlatIP(self.embedding_dim)
            self.documents = []
            self.doc_metadata = {}
            self.is_dirty = True
            logger.info("Cleared vector store")
        except Exception as e:
            logger.error(f"Error clearing vector store: {e}")

class VectorStoreManager:
    """Manages vector store instances per session"""
    
    def __init__(self):
        """Initialize vector store manager without S3 dependency"""
        self.sessions: Dict[str, FAISSVectorStore] = {}
        logger.info("Vector store manager initialized independently of S3")
    
    def get_or_create_store(self, session_id: str) -> Optional[FAISSVectorStore]:
        """Get or create a vector store for a session"""
        try:
            if not FAISS_AVAILABLE:
                logger.warning("FAISS not available, cannot create vector store")
                return None
            
            if session_id not in self.sessions:
                self.sessions[session_id] = FAISSVectorStore(session_id)
                logger.info(f"Created new vector store for session: {session_id}")
            
            return self.sessions[session_id]
            
        except Exception as e:
            logger.error(f"Error creating vector store for session {session_id}: {e}")
            return None
    
    # S3-specific method removed - use add_document_from_content instead
    
    def add_document_from_content(self, session_id: str, content: str, source_filename: str, metadata: dict = None) -> bool:
        """Add a document to vector store from content string - SIMPLE approach"""
        try:
            # Get vector store
            store = self.get_or_create_store(session_id)
            if not store:
                logger.error(f"❌ Failed to get/create vector store for session {session_id}")
                return False
            
            # Add to vector store directly
            chunks_added = store.add_document(content, source_filename, metadata or {})
            logger.info(f"✅ Added document {source_filename} to vector store: {chunks_added} chunks")
            return chunks_added > 0
            
        except Exception as e:
            logger.error(f"❌ Error adding document {source_filename} to vector store: {e}")
            import traceback
            logger.error(f"Full traceback: {traceback.format_exc()}")
            return False
    
    def add_document_from_file(self, session_id: str, file_path: str, source_filename: str) -> bool:
        """Add a document to vector store from local file"""
        try:
            # Get vector store
            store = self.get_or_create_store(session_id)
            if not store:
                return False
            
            # Read document from file
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            # Add to vector store
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
                # Use the new comprehensive method with much higher limits
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
    """Initialize the global vector store manager - S3 independent"""
    global vector_store_manager
    
    try:
        if not FAISS_AVAILABLE:
            logger.warning("FAISS not available, vector store disabled")
            return None
        
        vector_store_manager = VectorStoreManager()
        logger.info("Initialized vector store manager independently")
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
    return FAISS_AVAILABLE and vector_store_manager is not None