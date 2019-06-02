#!/bin/bash
if git grep --include "*.ipynb" '"image/png":' -R .; then
    >&2 echo "You may not check in binary data into the repository."
    exit 1
fi
