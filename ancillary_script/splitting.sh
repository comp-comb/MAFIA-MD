#!/bin/bash

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

[ ! -d "./split" ] && mkdir split

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
rm info lines times
