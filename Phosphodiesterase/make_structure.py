pdb = '/scratch/pyContacts/mdentropy/pde6/frame0.pdb'
out0 = 'dmi_side_output_not_norm_full_PDE6_list_0_stride10.pkl'
out_name = 'dmi_side_pde6_notnorm_avg_stride10_structure_over5.py'

import numpy as np
import mdtraj as md
import pickle


with open(out0, 'rb') as handle:
        out00 = pickle.load(handle)

out000 = np.copy(out00)


out000 -= out000.diagonal() * np.eye(*out000.shape)

t = md.load(pdb)
top = t.topology

a = np.where(out000 > 0.5)

first_frame = t.xyz

f = open(out_name,'w')
f.write('from pymol import cmd\n')
f.write('from pymol.cgo import *\n')
f.write("cmd.load('frame0.pdb', 'prot')\n")
f.write("cmd.show('cartoon')\n")
f.write('obj=[]\n')

for i in range(len(a[0])):
    j = a[0][i]
    k = a[1][i]
    l = top.select('resi ' + str(j) + ' and name CA')
    m = top.select('resi ' + str(k) + ' and name CA')
    n = first_frame[0,l,:] * 10
    o = first_frame[0,m,:] * 10
    c = out000[a[0][i],a[1][i]]
    c2 = np.round_(c, decimals=4)
    f.write('obj.extend([CYLINDER, ' + str(n[0][0]) + ', ' + str(n[0][1]) + ', ' + str(n[0][2]) + ', ' + str(o[0][0]) + ', ' + str(o[0][1]) + ', ' + str(o[0][2]) + ', ' + '0.15, ' + str(c2) + ', 0.0000, 0.0000, ' + str(c2) + ', 0.0000, 0.0000, ])\n' )

f.write("cmd.load_cgo(obj, 'contacts')")
f.close()
