#! /usr/bin/env bash

set -Eeuo pipefail

accessions_list_file="${1:-accessions.list}"
adapter_seq="TGGAATTCTCGGGTGCCAAGG"

timestamp="$(date +%s)"
runner="bin/final-thesis-${timestamp}"

function cleanup() {
    rm -rf "$runner"
}

trap cleanup ERR SIGINT SIGTERM SIGHUP SIGQUIT

CGO_ENABLED=0 go build -ldflags='-s -w -extldflags "-static"' -o "$runner"


while IFS= read -r accession_line; do
    accession=$(echo "${accession_line%#*}" | tr -d '[:space:]') # remove comments and trim spaces
    if [[ -z "$accession" ]]; then
        echo "empty accession line: ${accession_line}"
        continue
    fi
    "$runner" --accession "$accession" --adapter-seq "$adapter_seq"
done < "$accessions_list_file"
