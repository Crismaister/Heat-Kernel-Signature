#!/bin/false

function git-short-commit {
    git log -1 | head -n 1 | cut -d\  -f 2 | cut -b -8
}

function die {
    echo "$@"
    exit 1
}