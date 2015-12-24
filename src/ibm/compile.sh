#!/usr/bin/env bash

#array_nf=(1 2 3 4 5 10)
#array_nm=(1 2 3 4 5 10)
array_nf=(1 2 4)
array_nm=(1 2 4)

for i in "${array_nf[@]}"
do
    for j in "${array_nm[@]}"
    do
        echo $i;
        echo $j;
        myvar=$(cat age_vs_size_imprint.cpp | sed 's/const size_t Nfp = [0-9]/const size_t Nfp = '$i'/g')
        echo "$myvar" | sed 's/const size_t Nmp = [0-9]/const size_t Nmp = '$j'/g' > age_vs_size_imprint_nf"$i"_nm"$j".cpp
        g++ -Wall -ggdb -std=c++11 -o xage_vs_size_imprint_nf"$i"_nm"$j" age_vs_size_imprint_nf"$i"_nm"$j".cpp Individual.cpp -lm -lgsl -lgslcblas
        rm age_vs_size_imprint_nf"$i"_nm"$j".cpp
    done
done
