#!/bin/bash

rm math.out
#/Users/mjg/physics/research/software/./form 'boson_reduction.frm' >> math.out
~/physics/utils/form/./tform 'boson_reduction.frm' >> math.out
cat bos.tmp  | tr -d ' ' | tr -d '\n'
