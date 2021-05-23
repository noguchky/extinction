#!/bin/bash

SOURCEDIR=$(cd $(dirname $0) && pwd)

function show_usage() {
    echo "Usage' ./haddCoin.sh [emlist]"
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
coin_filenames=()

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
            coin_filename=${marged_filename/.root/_coin.root}
            # echo "marged_dirname = ${marged_dirname}"
            # echo "coin_filename  = ${coin_filename}"

            coin_filenames+=(${marged_dirname}/${coin_filename})
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
    coin_filename=${marged_filename/.root/_coin.root}
    # echo "marged_dirname = ${marged_dirname}"
    # echo "hists_filename = ${hists_filename}"
    # echo "spill_filename = ${spill_filename}"

    coin_filenames+=(${marged_dirname}/${coin_filename})
fi

echo "coin_filenames = ${coin_filenames[@]}"
if [ -z "${coin_filenames}" ]; then
    exit
fi

hadd_dirname=$(dirname ${coin_filenames[0]/marged/hadd})
hadd_coin_filename=fct_hadd_$(echo ${list} | grep -o -P '(?<=em_).+(?=\.txt)')_hists.root
echo ${hadd_dirname}
echo ${hadd_coin_filename}

if [ ! -d ${hadd_dirname} ]; then
    mkdir -p ${hadd_dirname}
fi

hadd -f ${hadd_dirname}/${hadd_coin_filename} ${coin_filenames[@]}

root -b -q ${SOURCEDIR}/haddCoin.C'("'${hadd_dirname}/${hadd_coin_filename}'")'
