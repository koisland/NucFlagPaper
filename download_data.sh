#!/bin/bash

set -euo pipefail

url=$1
output_dir=$2

url_bname="${url}"
output_dir_bname="${output_dir}"

wget --no-check-certificate "${url}" -P "${output_dir}"

if [[ ${output_dir_bname} == "assembly" ]] && [[ ${url_bname} == "*.gz" ]]; then
    gunzip "${output_dir}/${url_bname}"
fi
