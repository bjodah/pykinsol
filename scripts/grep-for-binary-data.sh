#!/bin/bash
if git grep '"image/png":' -R . --include "*.ipynb"; then
    >&2 echo "You may not check in binary data into the repository."
    exit 1
fi
