#!/bin/bash

places=(
    mito-10min_s1.php-40-47629
    hitachi-10min_a1.php-40-1011
)

dates=(
    # 2020-05-09
    # 2020-05-10
    # 2020-05-11
    # 2020-05-12
    # 2020-05-13
    # 2020-05-14
    # 2020-05-15
    # 2020-05-16
    # 2020-05-17
    # 2020-05-18
    # 2020-05-19
    # 2020-05-20
    # 2020-05-21
    # 2020-05-22
    # 2020-05-23
    # 2020-05-24
    # 2020-05-25
    # 2020-05-26
    # 2020-05-27
    # 2020-05-28
    # 2020-05-29
    # 2020-05-30
    # 2020-05-31

    # 2021-03-10
    # 2021-03-11
    # 2021-03-12
    # 2021-03-13
    # 2021-03-14
    # 2021-03-15
    # 2021-03-16
    # 2021-03-17
    # 2021-03-18
    # 2021-03-19
    # 2021-03-20
    # 2021-03-21
    # 2021-03-22
    # 2021-03-23
    # 2021-03-24
    # 2021-03-25
    # 2021-03-26
    # 2021-03-27
    # 2021-03-28
    # 2021-03-29
    # 2021-03-30
    # 2021-03-31
    # 2021-04-01
    # 2021-04-02
    # 2021-04-03
    # 2021-04-04
    # 2021-04-05
    2021-04-06
    2021-04-07
    2021-04-08
    2021-04-09
    2021-04-10
    2021-04-11
    2021-04-12
    2021-04-13
    2021-04-14
)

if ! [ -d ./data_jma ]; then
    mkdir ./data_jma
fi

for str1 in ${places[@]}; do
    vals=($(echo ${str1} | tr "-" " "))
    name=${vals[0]}
    page=${vals[1]}
    prefecture=${vals[2]}
    block=${vals[3]}

    for str2 in ${dates[@]}; do
        vals=($(echo ${str2} | tr "-" " "))
        year=${vals[0]}
        month=${vals[1]}
        day=${vals[2]}

        filename="${page}?prec_no=${prefecture}&block_no=${block}&year=${year}&month=${month}&day=${day}&view=p1"

        wget "http://www.data.jma.go.jp/obd/stats/etrn/view/${filename}" -O - \
            | grep 'tr class="mtx"' \
            | sed -e 's/ //g' \
            | sed -e 's/<[^>]*>/ /g' \
                  > ./data_jma/jma_${name}_${str2}.txt

        sleep 0.2
    done
done