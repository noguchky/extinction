#!/bin/bash

SOURCEDIR=$(cd $(dirname $0) && pwd)

function show_usage() {
    echo "Usage' ./haddHist.sh [emlist]"
    echo
}

show_events_option=
exec_single=0
while :; do
    if [ ${1}_ = "-h_" ]; then
        show_usage
        exit
    else
        break
    fi
done
list=${1}

if [ -z ${list} ]; then
    show_usage
    exit
elif [ ! -f "${list}" ]; then
    echo "emlist is not found"
    exit
fi

emcount=-1
boards=()
filenames=()
hists_filenames=()
spill_filenames=()

while read line; do
    echo "${line}"
    if [ -z "${line}" ]; then
        continue
    fi

    array=($(echo ${line}))
    this_emcount=${array[0]}
    this_board=${array[1]}
    this_date=${array[2]}
    this_filename=${array[3]}
    # echo "this_emcount  = ${this_emcount}"
    # echo "this_board    = ${this_board}"
    # echo "this_date     = ${this_date}"
    # echo "this_filename = ${this_filename}"

    this_root_dirname=$(dirname ${this_filename/data/ana})
    this_root_filename=$(basename ${this_filename/.dat/.root})
    # echo "this_root_dirname  = ${this_root_dirname}"
    # echo "this_root_filename = ${this_root_filename}"

    if [ ${emcount} -ne ${this_emcount} ]; then
        if [ ${#filenames[@]} -ne 0 ]; then
            marged_dirname=$(dirname ${filenames[0]//id[0-9][0-9][0-9][0-9]/marged})
            marged_filename=$(basename ${filenames[0]//id[0-9][0-9][0-9][0-9]/marged})
            hists_filename=${marged_filename/.root/_hists.root}
            spill_filename=${marged_filename/.root/_spill.root}
            # echo "marged_dirname  = ${output_dirname}"
            # echo "marged_filename = ${output_filename}"

            hists_filenames+=(${marged_dirname}/${hists_filename})
            spill_filenames+=(${marged_dirname}/${spill_filename})
        fi
        boards=()
        filenames=()
    fi

    emcount=${this_emcount}
    boards+=(${this_board})
    filenames+=(${this_root_dirname}/${this_root_filename})
done < <(cat ${list} | sort -k 1n,1 -k 2n,2)


if [ ${#filenames[@]} -ne 0 ]; then
    marged_dirname=$(dirname ${filenames[0]//id[0-9][0-9][0-9][0-9]/marged})
    marged_filename=$(basename ${filenames[0]//id[0-9][0-9][0-9][0-9]/marged})
    hists_filename=${marged_filename/.root/_hists.root}
    spill_filename=${marged_filename/.root/_spill.root}
    # echo "marged_dirname  = ${marged_dirname}"
    # echo "hists_filename = ${hists_filename}"
    # echo "spill_filename = ${spill_filename}"

    hists_filenames+=(${marged_dirname}/${hists_filename})
    spill_filenames+=(${marged_dirname}/${spill_filename})
fi

echo "hists_filenames = ${hists_filenames[@]}"
echo "spill_filenames = ${spill_filenames[@]}"
if [ -z "${hists_filenames}" ]; then
    exit
fi

hadd_dirname=$(dirname ${hists_filenames[0]/marged/hadd})
hadd_hists_filename=fct_hadd_$(echo ${list} | grep -o -P '(?<=em_).+(?=\.txt)')_hists.root
hadd_spill_filename=fct_hadd_$(echo ${list} | grep -o -P '(?<=em_).+(?=\.txt)')_spill.root
echo ${hadd_dirname}
echo ${hadd_hists_filename}
echo ${hadd_spill_filename}

if [ ! -d ${hadd_dirname} ]; then
    mkdir -p ${hadd_dirname}
fi

hadd -f ${hadd_dirname}/${hadd_hists_filename} ${hists_filenames[@]}
hadd -f ${hadd_dirname}/${hadd_spill_filename} ${spill_filenames[@]}

${SOURCEDIR}/../build/haddHist ${hadd_dirname}/${hadd_hists_filename}
