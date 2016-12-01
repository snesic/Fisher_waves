#!/bin/bash
#
# Run the 2d sFKPP equation in parallel.
# Creates a dir, named by the realization number, copies the program into the dir and runs a single iteration
#
clear

#system_size=128
#i_0=$2
#i_end=$3
#dir="$dir/$i_0"

# Input System size, start and end realization
echo "System size: "
read system_size
echo "Introduce the starting realization: "
read i_0
echo "Introduce the end realization: "
read i_end


for (( i = $i_0; i <= $i_end; i++ ))
do
#	echo $i
#	if [! -d "$i"]; then 
		mkdir $i
		echo "Dir $i created"
#	fi
    make
	cp fisher_waves $i/
    cd $i
	nohup ./fisher_waves $system_size > fisher_waves.out > fisher_waves.err < /dev/null &
	sleep 3
    cd ..
done


exit 0
