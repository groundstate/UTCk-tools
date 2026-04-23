#!/bin/bash
for file in original/*; do
	if [ -f "$file" ]; then
		echo "Processing $file"
		fbase=$(basename $file)
		echo $fbase
		editrnxnav.py --extract-header $file | grep -v IONOSPHERI > $fbase
	fi
done
