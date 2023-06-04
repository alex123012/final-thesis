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
    echo "running: ${accession}"
    "$runner" --accession "$accession" --adapter-seq "$adapter_seq"
done < "$accessions_list_file"

while IFS= read -r project; do
    project_samples=$(cat "$accessions_list_file" | grep "$project" | awk '{print $1}' | xargs | tr ' ' '|')
    files=$(ls results/tables/*.tsv | grep -E "$project_samples")

    echo "running: for project ${project}"
    python3 scripts/full_plots.py " of the $project project" "results/result-barplot-${project}.png" "results/result-heatmap-${project}.png" "results/result-histogram-${project}.png" "results/result-seed-eding-${project}.png" $files
done < <(awk '{print $3}' accessions.list | awk NF | sort -u)

echo "running for all samples"
python3 scripts/full_plots.py "" results/result-barplot-full.png results/result-heatmap-full.png results/result-histogram-full.png results/result-seed-eding-full.png results/tables/*.tsv
