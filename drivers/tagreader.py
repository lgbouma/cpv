"""
[WIP] - give some directory, get the Mac OS X Finder "Tags" associated with all
files in that directory.

Under the hood, this works as follows:

get_macosx_tags.sh:

    1) Pull binary extended attribute information using xattr
    2) Convert to an XML format
    3) Pull the "property list" (plist) information from the separate Mac-specific
        path that it lives at.
    4) Wrangle the resulting human-readable string names from the xmls.

This script is then the higher order wrapper.

A potential challenge right now is that this is ***slow*** - like 0.5 sec per file.
->Not great with thousands of files.  Might be multithreadable?
"""
import numpy as np, pandas as pd
import os
from glob import glob
import subprocess

def get_file_tags(file_path):
    script_path = "./get_macosx_tags.sh"
    output = subprocess.check_output([script_path, file_path],
                                     stderr=subprocess.STDOUT,
                                     universal_newlines=True)
    tags = output.strip().split("\n")
    if 'No such xattr' in tags[0]:
        tags = ''
    return tags

def get_directory_tags(directory):
    file_tags = {}
    filepaths = np.sort(glob(os.path.join(directory, '*pdf')))

    for ix, file_path in enumerate(filepaths):
        tags = get_file_tags(file_path)
        print(ix, file_path, tags)
        if tags:
            # Remove numeric tags and convert remaining tags to a list
            tags = [tag for tag in tags if not tag.isdigit()]
            file_tags[file_path] = tags
    return file_tags

# Example usage
directory_path = "/Users/luke/vet_subset"
tags_by_file = get_directory_tags(directory_path)

# Create a DataFrame from the file paths and tags
df = pd.DataFrame(list(tags_by_file.items()), columns=['File', 'Tags'])

# Print the DataFrame
df.to_csv('vet_subset.csv', index=False)
