import argparse
import os
from pathlib import Path
from langchain_community.document_loaders import PyPDFLoader
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.vectorstores import FAISS
from langchain_aws import BedrockEmbeddings
from langchain_aws import ChatBedrock
from langchain.chains import create_retrieval_chain
from langchain.chains.combine_documents import create_stuff_documents_chain
from langchain_core.prompts import ChatPromptTemplate

def process_pdfs(folder_path):
    """Process PDF files in a folder and create/load FAISS vector store"""
    folder_name = Path(folder_path).name
    index_name = f"{folder_name}.faiss"
    index_path = Path(folder_path) / index_name
    
    pdf_files = [f for f in os.listdir(folder_path) if f.lower().endswith('.pdf')]
    
    if not pdf_files:
        raise ValueError("No PDF files found in the specified folder")
    
    # Check for existing index
    if index_path.exists():
        print(f"Loading existing FAISS index: {index_path}")
        return FAISS.load_local(folder_path, BedrockEmbeddings(), index_name=folder_name)
    
    print(f"Processing {len(pdf_files)} PDF files...")
    
    # Process all PDFs
    text_splitter = RecursiveCharacterTextSplitter(
        chunk_size=500,
        chunk_overlap=100
    )
    
    all_chunks = []
    for pdf_file in pdf_files:
        loader = PyPDFLoader(str(Path(folder_path) / pdf_file))
        pages = loader.load()
        chunks = text_splitter.split_documents(pages)
        all_chunks.extend(chunks)
    
    # Create and save vector store
    embeddings = BedrockEmbeddings()
    vector_store = FAISS.from_documents(all_chunks, embeddings)
    vector_store.save_local(folder_path, index_name=folder_name)
    print(f"Created FAISS index at: {index_path}")
    return vector_store

def main():
    parser = argparse.ArgumentParser(description='PDF Question Answering System')
    parser.add_argument('folder', type=str, help='Path to folder containing PDF files')
    args = parser.parse_args()

    # Process PDFs and get vector store
    vector_store = process_pdfs(args.folder)
    retriever = vector_store.as_retriever()
    
    # Set up Bedrock LLM with Claude 3
    llm = ChatBedrock(
        model_id="anthropic.claude-3-sonnet-20240229-v1:0",
        model_kwargs={
            "temperature": 0.5,
            "max_tokens": 2048
        },
        region_name="us-west-2"
    )
    
    # Modified prompt template for Claude
    prompt = ChatPromptTemplate.from_template(
        """\n\nHuman: You are a helpful assistant. Answer the question based only on this context:
        <context>
        {context}
        </context>
        
        Question: {input}
        
        Assistant:"""
    )
    question_answer_chain = create_stuff_documents_chain(llm, prompt)
    qa_chain = create_retrieval_chain(retriever, question_answer_chain)
    
    # Start interactive Q&A
    print("\nEnter questions about the PDF content (type 'exit' to quit):")
    while True:
        query = input("\nQuestion: ")
        if query.lower() in ('exit', 'quit'):
            break
        response = qa_chain.invoke({"input": query})
        print(f"\nAnswer: {response['answer']}")

if __name__ == "__main__":
    main()
