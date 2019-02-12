import numpy as np
import pickle
import matplotlib.pyplot as plt
import msmexplorer as msme

obj = 'dmi_side_PDE6.pkl'
out = 'dmi_side_output_PDE6.pkl'
out_name = 'dmi_side_pde6'
o_type = 'pdf'

with open(obj, 'rb') as handle:
        obj1 = pickle.load(handle)

with open(out, 'rb') as handle:
        out00 = pickle.load(handle)


out000 = np.copy(out00)

out000 -= out000.diagonal() * np.eye(*out000.shape)

lab_arr = obj1.labels + 1
lab_list = lab_arr.tolist()
blank = [''] * len(lab_list)
blank[2::20] = lab_list[2::20]

ax = msme.plot_chord(out000, labels=blank, threshold=0.4)
plt.savefig(out_name + '_chord.' + o_type)
plt.clf()

#out1 -= out1.diagonal() * np.eye(*out1.shape)
out000[out000 < 0.05] = np.nan

minr = lab_arr.min()
maxr = lab_arr.max()

plt.imshow(out000, origin='lower', interpolation='none', extent=[minr, maxr, minr, maxr], cmap='viridis_r')
plt.colorbar()
plt.xlabel('Residue Number', fontsize=16)
plt.ylabel('Residue Number', fontsize=16)
plt.savefig(out_name + '_matrix.' + o_type)
plt.clf()

