#!/usr/bin/env python3

import numpy as np
import datetime

# how rapidly mortality accelerates with increased care
c = [ 0.01]

# whether slope of kin competition function is concave, convex or linear  
a = [ 0.5]

maxs = [100]
 
td = [ .1 ]

se = [ 3 ]

mmin = [ 0.02 ] 

u = [ 0.5, 0.5 ]

v = [ 1, 1]

lambda_val = 1.0

b = [ 0.5, 0.3]

r_mother_brood = 0.5
r_mother_off = 0.5

ctr = 0

exe = "care_transgen.exe"

def space_str(vec):
    vect = [ str(x) for x in vec ]

    vecs = " ".join(vect)

    return(vecs)

date = datetime.datetime.now()
base_name = "iter_care_transgen_" +\
        f"{date:%d}_{date:%m}_{date:%Y}_{date:%H}{date:%M}{date:%S}"

for c_i in c:
    for a_i in a:
        for maxs_i in maxs:
            for td_i in td:
                for se_i in se:
                    for mmin_i in mmin:

                        m = [ mmin_i * 2, mmin_i * 2, mmin_i * 2]

                        print("echo " + str(ctr))

                        print("./" + exe + " " \
                                + str(c_i) + " " \
                                + str(a_i) + " " \
                                + str(maxs_i) + " " \
                                + str(td_i) + " " \
                                + str(se_i) + " " \
                                + str(mmin_i) + " " \
                                + space_str(u) + " " \
                                + space_str(v) + " " \
                                + str(lambda_val) + " " \
                                + space_str(b) + " " \
                                + space_str(m) + " " \
                                + str(r_mother_brood) + " " \
                                + str(r_mother_off) + " " \
                                + base_name + "_" + str(ctr))

                        ctr+=1






