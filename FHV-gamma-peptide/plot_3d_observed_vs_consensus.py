"""
plot 3d pmf
plot observed and consensus distrabution for highest divergence window
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt

active =  pickle.load(open('active_set_multi_e2_full_mix.pickle', 'rb'))
con1 = pickle.load(open('div_c_78_multi_e4_full123_mix.pickle', 'rb'))
obs1 = pickle.load(open('div_o_78_multi_e4_full123_mix.pickle', 'rb'))

tic_cl = pickle.load(open('tic_cl_full_centers3_try2.pickle', 'rb'))
final = pickle.load(open('tram_mix_3dim_f_full_state_mix.pickle', 'rb'))
cl = pickle.load(open('clust_col_skip_obj_cl_full.pickle', 'rb'))

f_kcal = final * 0.6

reshape = np.reshape(f_kcal, (10, 451), order='F')

p = cl.clustercenters / 16
x, y = np.meshgrid(p, tic_cl[:,0])
x2, y2 = np.meshgrid(p, tic_cl[:,1])

# 3d pmf CV vs IC2
fig1, ax1 = plt.subplots()
cs = ax1.scatter(x, y, s=100, c=reshape, vmax=19.5)
ax1.set_xlabel('Percent Helicity', fontsize=16)
ax1.set_ylabel('IC 2', fontsize=16)
cbar = fig1.colorbar(cs)
cbar.set_label('PMF (kcal/mol)', fontsize=16)
fig1.tight_layout()
fig1.savefig('2d_ic2_mix_dots.eps')

# 3d pmf CV vs IC3
fig2, ax2 = plt.subplots()
cs = ax2.scatter(x2, y2, s=100, c=reshape, vmax = 19.5)
ax2.set_xlabel('Percent Helicity', fontsize=16)
ax2.set_ylabel('IC 3', fontsize=16)
cbar = fig2.colorbar(cs)
cbar.set_label('PMF (kcal/mol)', fontsize=16)
fig2.tight_layout()
fig2.savefig('2d_percentv3_mix_dots.eps')

zer_o = np.zeros(f_kcal.shape)
zer_c = np.zeros(f_kcal.shape)
zer_o[active] = obs1
zer_c[active] = con1
re_o = np.reshape(zer_o, (10, 451), order='F')
re_c = np.reshape(zer_c, (10, 451), order='F')
re_o[re_o < 0.0001] = np.nan
re_c [re_c < 0.0001] = np.nan

# 3d Observed distribution of max divergence window CV vs IC2
fig3, ax3 = plt.subplots()
cs = ax3.scatter(x, y, s=100, c=re_o, vmax=0.143)
ax3.set_xlabel('Percent Helicity', fontsize=16)
ax3.set_ylabel('IC 2', fontsize=16)
ax3.set_title('Win 72 Observed distribution', fontsize=16)
ax3.set_ylim(-1.5,0.3)
ax3.set_xlim(0.47,0.52)
cbar = fig3.colorbar(cs)
cbar.set_label('Probability', fontsize=16)
fig3.tight_layout()
fig3.savefig('2d_ic2_mix_dots_observed.eps')

#3d Consensus distribution of max divergence window CV vs IC2
fig4, ax4 = plt.subplots()
cs = ax4.scatter(x, y, s=100, c=re_c, vmax=0.143)
ax4.set_xlabel('Percent Helicity', fontsize=16)
ax4.set_ylabel('IC 2', fontsize=16)
ax4.set_title('Win 72 Consensus distribution', fontsize=16)
ax4.set_ylim(-1.5,0.3)
ax4.set_xlim(0.47,0.52)
cbar = fig4.colorbar(cs)
cbar.set_label('Probability', fontsize=16)
fig4.tight_layout()
fig4.savefig('2d_ic2_mix_dots_consensus.eps')

#3d Observed distribution of max divergence window CV vs IC3
fig5, ax5 = plt.subplots()
cs = ax5.scatter(x2, y2, s=100, c=re_o2, vmax=0.143)
ax5.set_xlabel('Percent Helicity', fontsize=16)
ax5.set_ylabel('IC 3', fontsize=16)
ax5.set_title('Win 72 Observed distribution', fontsize=16)
ax5.set_ylim(-1.05,2.4)
ax5.set_xlim(0.47,0.52)
cbar = fig5.colorbar(cs)
cbar.set_label('Probability', fontsize=16)
fig5.tight_layout()
fig5.savefig('2d_ic3_mix_dots_observed.eps')

#3d Consensus distribution of max divergence window CV vs IC3
fig6, ax6 = plt.subplots()
cs = ax6.scatter(x2, y2, s=100, c=re_c2, vmax=0.143)
ax6.set_xlabel('Percent Helicity', fontsize=16)
ax6.set_ylabel('IC 3', fontsize=16)
ax6.set_title('Win 72 Consensus distribution', fontsize=16)
ax6.set_ylim(-1.05,2.4)
ax6.set_xlim(0.47,0.52)
cbar = fig6.colorbar(cs)
cbar.set_label('Probability', fontsize=16)
fig6.tight_layout()
fig6.savefig('2d_ic3_mix_dots_consensus.eps')
