#!/bin/sh

ME=`basename "$0"`
if /sbin/pidof -o %PPID -x "$ME">/dev/null; then
    echo "Script $ME already running, skipping"
    exit 0
fi


cd "$(dirname "$0")"/../..
python ./scripts/cleanup_datasets/cleanup_datasets.py -d 60 -5 -r "$@" >> ./scripts/cleanup_datasets/purge_folders.log
