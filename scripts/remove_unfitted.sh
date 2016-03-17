#!/bin/bash

verbose=true
if [[ ( $1 == "--silent") ||  $1 == "-s" ]] 
then 
  verbose=false
fi

print_files() { 
	echo -e "\nThe following files will be removed:" 
  arg=("$@")

  for ((i=0;i<=$#;i++)) ; do
    echo "${arg[i]}"
  done
	} 

suffix="_BEAGLE"

list="$(find . -name '*fits.gz' -size 0)"

for file in $list ; do
  i=$(awk -v a="$file" -v b="${suffix}.fits.gz" 'BEGIN{print index(a,b)}')

  file=${file:0:$i-1}

  rm_files="${file}${suffix}.fits.gz ${file}.lock ${file}${suffix}_MN*"

  if [ "$verbose" = true ] ; then
    print_files ${rm_files}
  fi

  rm ${rm_files}
done  
