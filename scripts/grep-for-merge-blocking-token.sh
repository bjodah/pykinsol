#!/bin/bash
if git grep "DO-NOT-MERGE!" -R . --exclude $(basename $0); then
    exit 1
fi
