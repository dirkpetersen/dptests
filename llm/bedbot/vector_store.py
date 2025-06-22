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
                 chunk_size: int = 512, chunk_overlap: int = 50):
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
        """Split text into overlapping chunks"""
        chunks = []
        
        # Simple sentence-aware chunking
        sentences = text.split('. ')
        current_chunk = ""
        chunk_count = 0
        
        for sentence in sentences:
            # Add sentence to current chunk
            test_chunk = current_chunk + sentence + ". "
            
            if len(test_chunk) <= self.chunk_size:
                current_chunk = test_chunk
            else:
                # Current chunk is full, save it and start new one
                if current_chunk.strip():
                    chunk_id = f"{source}_{chunk_count}"
                    chunks.append(DocumentChunk(
                        text=current_chunk.strip(),
                        source=source,
                        chunk_id=chunk_id,
                        metadata={'chunk_index': chunk_count}
                    ))
                    chunk_count += 1
                
                # Start new chunk with overlap
                if self.chunk_overlap > 0 and current_chunk:
                    # Take last chunk_overlap characters as overlap
                    overlap_text = current_chunk[-self.chunk_overlap:]
                    current_chunk = overlap_text + sentence + ". "
                else:
                    current_chunk = sentence + ". "
        
        # Add final chunk if not empty
        if current_chunk.strip():
            chunk_id = f"{source}_{chunk_count}"
            chunks.append(DocumentChunk(
                text=current_chunk.strip(),
                source=source,
                chunk_id=chunk_id,
                metadata={'chunk_index': chunk_count}
            ))
        
        logger.info(f"Chunked document {source} into {len(chunks)} chunks")
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
            # Generate document hash for deduplication
            doc_hash = hashlib.md5(f"{source}_{text}".encode()).hexdigest()
            
            # Check if document already exists
            if doc_hash in self.doc_metadata:
                logger.info(f"Document {source} already exists in vector store")
                return 0
            
            # Chunk the document
            chunks = self._chunk_text(text, source)
            
            if not chunks:
                logger.warning(f"No chunks generated for document: {source}")
                return 0
            
            # Generate embeddings for all chunks
            chunk_texts = [chunk.text for chunk in chunks]
            embeddings = self.embedding_model.encode(chunk_texts, convert_to_numpy=True)
            
            # Normalize embeddings for cosine similarity
            faiss.normalize_L2(embeddings)
            
            # Add to FAISS index
            self.index.add(embeddings)
            
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
                'timestamp': datetime.now().isoformat()
            }
            
            self.is_dirty = True
            logger.info(f"Added document {source} with {len(chunks)} chunks to vector store")
            return len(chunks)
            
        except Exception as e:
            logger.error(f"Error adding document {source} to vector store: {e}")
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
    
    def get_all_documents_context(self, max_chars: int = 15000) -> str:
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
            logger.info(f"Processing {num_documents} documents for comprehensive analysis")
            
            # Calculate chars per document to ensure all documents are included
            chars_per_doc = min(500, max_chars // num_documents) if num_documents > 0 else 500
            logger.info(f"Allocating {chars_per_doc} chars per document to fit all {num_documents} documents")
            
            # Take one chunk from each document to ensure broad coverage
            for source, chunks in by_source.items():
                # Take the first chunk from each document (usually contains name/header info)
                chunk = chunks[0]
                source_text = f"\n--- From {source} ---\n"
                
                # Limit chunk text to allocated chars per document
                chunk_text = chunk.text[:chars_per_doc]
                if len(chunk.text) > chars_per_doc:
                    chunk_text += "..."
                chunk_text += "\n"
                
                context_parts.append(source_text)
                context_parts.append(chunk_text)
                total_chars += len(source_text) + len(chunk_text)
            
            context = "".join(context_parts).strip()
            logger.info(f"Generated comprehensive context: {len(context)} chars from {len(by_source)} documents")
            return context
            
        except Exception as e:
            logger.error(f"Error generating comprehensive context: {e}")
            return ""

    def get_context_for_query(self, query: str, max_chunks: int = 5, max_chars: int = 4000) -> str:
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
        """Get vector store statistics"""
        active_docs = sum(1 for meta in self.doc_metadata.values() if not meta.get('removed', False))
        total_chunks = sum(meta['chunk_count'] for meta in self.doc_metadata.values() if not meta.get('removed', False))
        
        return {
            'session_id': self.session_id,
            'total_documents': active_docs,
            'total_chunks': total_chunks,
            'index_size': self.index.ntotal,
            'embedding_model': self.embedding_model_name,
            'embedding_dim': self.embedding_dim,
            'is_dirty': self.is_dirty
        }
    
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
    
    def __init__(self, s3_client=None, s3_bucket=None):
        self.sessions: Dict[str, FAISSVectorStore] = {}
        self.s3_client = s3_client  
        self.s3_bucket = s3_bucket
    
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
    
    def add_document_from_s3(self, session_id: str, s3_key: str, source_filename: str) -> bool:
        """Add a document to vector store from S3"""
        try:
            if not self.s3_client or not self.s3_bucket:
                logger.error("S3 client or bucket not configured")
                return False
            
            # Get vector store
            store = self.get_or_create_store(session_id)
            if not store:
                return False
            
            # Read document from S3
            response = self.s3_client.get_object(Bucket=self.s3_bucket, Key=s3_key)
            content = response['Body'].read().decode('utf-8', errors='ignore')
            
            # Add to vector store
            chunks_added = store.add_document(content, source_filename, {'s3_key': s3_key})
            logger.info(f"Added document {source_filename} from S3 to vector store: {chunks_added} chunks")
            return chunks_added > 0
            
        except Exception as e:
            logger.error(f"Error adding document from S3 to vector store: {e}")
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
    
    def get_all_documents_for_session(self, session_id: str, max_chars: int = 15000) -> str:
        """Get context from all documents in a session's vector store"""
        try:
            if session_id not in self.sessions:
                logger.info(f"No vector store found for session: {session_id}")
                return ""
            
            return self.sessions[session_id].get_all_documents_context(max_chars)
            
        except Exception as e:
            logger.error(f"Error getting all documents for session {session_id}: {e}")
            return ""

    def get_context_for_session(self, session_id: str, query: str, max_chunks: int = 5, max_chars: int = 4000) -> str:
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

# Global vector store manager instance
vector_store_manager = None

def initialize_vector_store_manager(s3_client=None, s3_bucket=None) -> Optional[VectorStoreManager]:
    """Initialize the global vector store manager"""
    global vector_store_manager
    
    try:
        if not FAISS_AVAILABLE:
            logger.warning("FAISS not available, vector store disabled")
            return None
        
        vector_store_manager = VectorStoreManager(s3_client, s3_bucket)
        logger.info("Initialized vector store manager")
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