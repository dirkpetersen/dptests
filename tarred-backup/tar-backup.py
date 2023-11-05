#! /usr/bin/env python3

import sys, os, tarfile, shutil, tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed

def get_file_size(filename):
    """Returns the size of the file in bytes."""
    if os.path.isfile(filename):
        return os.path.getsize(filename)
    return 0

def tar_and_compress(source_dir, temp_dir, max_size=2**30):  # 1 GiB in bytes
    """
    Tars and compresses the files in a directory using gzip compression,
    ensuring that the tarball does not exceed the specified max size.
    """
    part = 1
    tarball_path = os.path.join(temp_dir, f"{os.path.basename(source_dir)}_part{part}.tar.gz")
    tar = tarfile.open(tarball_path, "w:gz")
    current_size = 0
    
    for dirpath, _, filenames in os.walk(source_dir):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)
            file_size = get_file_size(file_path)
            
            # Check if adding the next file would exceed the max size limit
            if current_size + file_size > max_size:
                # Close the current tar file
                tar.close()
                # Copy the tar file to the original directory
                shutil.copy2(tarball_path, source_dir)
                # Start a new tar file
                part += 1
                tarball_path = os.path.join(temp_dir, f"{os.path.basename(source_dir)}_part{part}.tar.gz")
                tar = tarfile.open(tarball_path, "w:gz")
                current_size = 0  # Reset the size counter for the new tar file
            
            # Add the file to the tar file
            tar.add(file_path, arcname=os.path.relpath(file_path, start=source_dir))
            current_size += file_size
    
    # Close and copy the last tar file
    tar.close()
    shutil.copy2(tarball_path, source_dir)

def process_directory(temp_dir, directory, max_size):
    """
    Processes a single directory to create compressed tar.gz files.
    """
    tar_and_compress(directory, temp_dir, max_size)

def main(input_directory, max_parallel=10, max_size=2**30):  # 1 GiB in bytes
    with tempfile.TemporaryDirectory() as temp_dir:
        directories = [os.path.join(input_directory, d) for d in os.listdir(input_directory)
                       if os.path.isdir(os.path.join(input_directory, d))]

        with ThreadPoolExecutor(max_workers=max_parallel) as executor:
            futures = [executor.submit(process_directory, temp_dir, d, max_size) for d in directories]

            for future in as_completed(futures):
                try:
                    future.result()
                    print(f"Completed processing for directory.")
                except Exception as exc:
                    print(f"Generated an exception: {exc}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: script.py <input_directory>")
        sys.exit(1)
    
    input_directory = sys.argv[1]  # Take the first argument as input directory
    max_size = 2**30  # Maximum size in bytes for each tar.gz file
    main(input_directory, max_parallel=10, max_size=max_size)
