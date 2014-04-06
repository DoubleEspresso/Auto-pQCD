#!/bin/bash

rm dat.txt
~/physics/Research/software/./form tmp.frm >> dat.txt

cat dat.txt | sed '/FORM/,/\.end/d' | sed '$d' | sed '$d' 

#expr=$(/Users/mikeglatzmaier/Desktop/Physics/Research/CalcTools/FORM/./form tmp.frm | sed '/FORM/,/.end/d' | sed '$d' | sed '$d')

#echo ${expr%?}
