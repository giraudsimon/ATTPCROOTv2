#!/bin/sh
# MODE=0 Evt2Raw followed by Raw2Cal; MODE=1 Evt2Cal
EVT_DIR=/mnt/simulations/ceclub/giraud/attpc/ATTPCROOTv2/macro/Unpack_HDF5/hdf5Files/
ROOT_DIR=/mnt/simulations/ceclub/giraud/attpc/ATTPCROOTv2/macro/Unpack_HDF5/rootFiles/
BIN_DIR="."
file_list=($(ls -t ${EVT_DIR}))

while : ; do 
    	file_list_size=${#file_list[@]}
    	cur_files=($(ls -t ${EVT_DIR}))
	cur_files_size=${#cur_files[@]}
	let "size_diff = $cur_files_size - $file_list_size"
    	if [[ $size_diff -gt 0 ]]; then
		echo "added" $size_diff "file(s)"
           	file_list=("${cur_files[@]}")
		for ((i=0; i<$size_diff; i++)) ; do
			echo "${cur_files[$i]}"
			runNum=$(echo "${cur_files[$i]}" | egrep -o '[[:digit:]]{4}' | head -n1)
			root -q "unpack.C("$runNum")"
			echo " "
			done
	elif [[ $size_diff -lt 0 ]]; then
		echo "removed $size_diff file(s) : " | sed 's/-//' 
		echo ${file_list[@]} ${cur_files[@]} | sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g" | tr ' ' '\n' | sort | uniq -u
		echo " "
           	file_list=("${cur_files[@]}") 
    	fi

    #echo "Waiting for changes."
    sleep $(expr 10 \* 2)
done
