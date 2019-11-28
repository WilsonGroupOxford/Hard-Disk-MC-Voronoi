#!/bin/bash

# vmd file to generate
export vmdFile=colloid.vmd

# file containing xyz information
export xyzFile=test.xyz
export diaFile=test_dia.dat

# header for vmd file
echo "package require topotools 1.5" > $vmdFile
echo "mol new test.xyz" >> $vmdFile 

line=1
num=$(head -1 $xyzFile)
for i in `seq 0 1 $num`
do
    r=$(head -$line $diaFile | tail -1 | awk '{print$1}')
    echo $r
    printf -v selection "%s%s%s%s%s" "name Ar" $i
    echo "mol addrep 0" >> $vmdFile
    echo "mol modselect $i" 0 $selection >> $vmdFile
    echo "mol modstyle" $i 0 "CPK" $r "0 40 0" >> $vmdFile
    line=$[$line+1]
done
echo "topo clearbonds" >> $vmdFile

# launch vmd
vmd -startup colloid.vmd

#menu main on
