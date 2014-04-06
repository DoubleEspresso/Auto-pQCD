#!/usr/bin/python
import os
import subprocess



MTH_CMD="MathKernel -script"
MTH_FILE="~/physics/research/proton_spin/form_codes/v1/src/integration/lattice_integrate.m"



print os.system(MTH_CMD + " " + MTH_FILE)
