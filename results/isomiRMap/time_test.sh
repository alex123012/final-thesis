#! /usr/bin/env bash

set -Eeuo pipefail

ISOMIRMAP_PATH="${1}"
FASTQ_FILES_PATH="${2}"

mkdir -p result time
function run_time(){
    time=$(which time)
    $time -a -o "${1}" python3 "$ISOMIRMAP_PATH/IsoMiRmap.py" ${2} --m "$ISOMIRMAP_PATH/MappingBundles/miRBase" --p "${3}" "${4}"

}

function process_time_file(){
   grep -v 'pagefaults' "${1}" | awk '{print $1}' | sed 's/user$//g'
}

echo -e "vanilla\tediting\tsample" | tee result_time.tsv >/dev/null
for file in "$FASTQ_FILES_PATH"/*; do
    echo "running for $file file"
    result_prefix_tmp=$(basename "$file")
    result_prefix="${result_prefix_tmp%.*}"
    result_prefix="${result_prefix%.*}"

    mkdir -p "result/$result_prefix"
    # rm -rf "time/${result_prefix}*.txt"
    # for _ in {1..3}; do
    #     run_time "time/${result_prefix}_vanilla.txt" "" "result/$result_prefix/out" "$file"
    #     run_time "time/${result_prefix}_editing.txt" "--e c" "result/$result_prefix/out_edit" "$file"
    # done

    paste \
        <(process_time_file "time/${result_prefix}_vanilla.txt") \
        <(process_time_file "time/${result_prefix}_editing.txt") \
        <(printf -- "$result_prefix\n%.0s" {1..3}) \
        --delimiters '\t' \
      | tee -a result_time.tsv >/dev/null
done

mean_time_comparasion=$(python3 -c 'import pandas as pd; df = pd.read_csv("result_time.tsv", sep="\t").dropna().groupby("sample").mean(); print((df.editing / df.vanilla).mean())')
echo "Editing takes longer than vanilla: for: $mean_time_comparasion"
sed -i "s#<!-- placeholder -->.*<!-- end placeholder -->#<!-- placeholder -->\`${mean_time_comparasion}\`<!-- end placeholder -->#g" README.md
