
# Assume the filename is 'source_file.py' and is in the same directory.
filename = 'source_file.py'

# Initialize an empty dictionary where the toolchain info will be stored
toolchain = {}

with open(filename, 'r') as file:
    # Read the file content
    file_content = file.read()
    # Execute the content in the context of the 'toolchain' dictionary
    exec(file_content, {'toolchain': toolchain})

# At this point, 'toolchain' should have the contents of the dictionary
# that was in the source file, assuming that the file contains a line like:
# toolchain = {'name': 'intel', 'version': '2020b'}
print(toolchain)

