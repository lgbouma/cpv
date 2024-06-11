#!/bin/bash

#  Description:
#  Please see docstring to tagreader.py.

if [ $# -eq 0 ]; then
    echo "Please provide a file path as an argument."
    exit 1
fi

file_path="$1"

if [ ! -f "$file_path" ]; then
    echo "File not found: $file_path"
    exit 1
fi

tags=$(xattr -p com.apple.metadata:_kMDItemUserTags "$file_path" | xxd -p -r | plutil -convert xml1 - -o - | xpath -e '/plist/array/string/text()' 2>/dev/null)

echo "$tags"
