#!/bin/bash

# bash script to convert a large trajectory file (xyz) containing multiple timesteps into a number of small xyz file containing only one timestep.
# the output files will be stored in a separated directory named "split"

# use: ./splitting.sh <name_of_large_trajectory_file>

grep -n Atoms $1| cut -d : -f 1 >lines 
grep -n Atoms $1 | cut -d : -f 3 >times

paste lines times > info
lineNumber=()
Timestep=()


while IFS= read -r line; do
    lineNumber+=("$line")
done < lines

while IFS= read -r line; do
    Timestep+=("$line")
done < times


a=$(cat lines|wc -l)
b=$(cat $1 |wc -l)
for (( i=0; i<$a; i++ ))
do
    ((x= ${lineNumber[$i]} - 1 ))
    ((y= ${lineNumber[$i + 1]} - 2 ))
    ((c= a -1))
    ((d= b ))
    ((name= ${Timestep[$i]} ))

    if [[ $i != $c ]] 
    then 
        sed -n ' '$x' , '$y' p ' $1 >./split/split_time_${name}.xyz
    else
        sed -n ' '$x' , '$d' p ' $1 >./split/split_time_${name}.xyz
    fi
done
