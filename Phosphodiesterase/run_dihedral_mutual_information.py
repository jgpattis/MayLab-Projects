import mdentropy.metrics
import mdtraj as md
import numpy as np
import pickle

t = md.load('frame0.pdb')
top = t.top
trj = md.load('PDE6_APO_no_water.xtc', top=top, stride=10)

dmi = mdentropy.metrics.DihedralMutualInformation(types=['phi', 'psi', 'chi1', 'chi2', 'chi3', 'chi4'], n_bins=3, method='knn', normed=False)
out = dmi.partial_transform(trj)

print('entropy calulation is complete for dihedral plus sidechain')
print('output is ', out)

p = pickle.HIGHEST_PROTOCOL
with open('dmi_side_PDE6_.pkl','wb') as handle:
        pickle.dump(dmi, handle,protocol=p)

with open('dmi_side_output_PDE6.pkl','wb') as handle:
        pickle.dump(out, handle,protocol=p)

x = np.copy(out)
np.savetxt('six_dmi.dat', x, fmt='%12.4f')
