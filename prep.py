import argparse
import os

# Set up argument parser
parser = argparse.ArgumentParser(description="Process a directory path.")
parser.add_argument("directory", type=str, help="Path to the directory")

# Parse arguments
args = parser.parse_args()

# Validate the directory
if not os.path.isdir(args.directory):
    print(f"Error: {args.directory} is not a valid directory")
else:
    print(f"Processing directory: {args.directory}")
    os.chdir(args.directory)

    for i in os.listdir():
        if '.ipynb' in i and 'checkpoint' not in i:
            cmd = 'jupyter nbconvert --to script ' + i
            os.system(cmd) 