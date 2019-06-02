#!/bin/bash
if git grep --exclude $(basename $0) "DO-NOT-MERGE!" -R .; then
    exit 1
fi
