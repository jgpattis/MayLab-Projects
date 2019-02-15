import numpy as np
import pickle
import pyemma.plots as mplt
import matplotlib
import matplotlib.pyplot as plt
## set path
indir = '/scratch/shn13007/mixfolder'
tica_file = '/jason_tram/mixTICs-ca-con-8.pkl'
## outfile names
out = 'mix_ca_con_8_red'
ics = 3  # look at this many ICs
edge = False
## load colvars 
colvar_list = [indir + "/comboCOLVAR{:2.1f}".format(i) for i in np.arange(0.7,15.5,0.1)]
col = [np.loadtxt(f, skiprows=1) for f in colvar_list]
## load tica info
tica = pickle.load(open(indir + tica_file, 'rb'))
for i in tica:
	i[:,0] *= -1
col1 = [ i[:]/16 for i in col ]
zipped = tuple(zip(col1, tica))
zipped2 = tuple(zip(tica, col1))
zipped3 = tuple(zip(tica,col))
## define plotting

centers_path = '/scratch/shn13007/mixfolder/jason_tram/cluster/thirty_100.pkl'
cl = pickle.load(open(centers_path, 'rb'))
cl1 = cl.clustercenters

def set_pub():
	plt.rc('xtick', labelsize=24)
	plt.rc('ytick', labelsize=24)
	plt.rc('lines', lw=1, color='k')

def plotFES(data,ic,outfile):
	""" plot all data as FES
	data: TICA data
	ic: second ic to plot against
	outfile: name of plot
	"""
	set_pub()
	fig1 = plt.figure(figsize=(8,6))
	fig1, ax1 = mplt.plot_free_energy(np.vstack(data)[:,0], np.vstack(data)[:,ic], kT=0.596, cbar_label='Free energy (kcal/mol)');
	ax1.set_xlabel('tIC 1', fontsize=24)
	ax1.set_ylabel('tIC ' + str(ic + 1), fontsize=24)
	plt.scatter(cl1[20,0],cl1[20,ic], marker="$1$", s=100, c='r')
	plt.scatter(cl1[17,0],cl1[17,ic],marker="$2$", s=100, c='r')
	plt.scatter(cl1[3,0],cl1[3,ic],marker="$3$", s=100, c='r')
	plt.scatter(cl1[18,0],cl1[18,ic],marker="$4$", s=100, c='r')
	plt.tight_layout()
	fig1.savefig(outfile + 'fes_1v' + str(ic + 1) + '.eps')
	fig1.clf()

def checklin(z,ic,outfile):
	""" plot US CV vs tIC1 to check they are linear
	z: give zipped colvar file and tica
	ic: give the tIC you want to plot
	outfile: name of plot
	"""
	set_pub()
	fig = plt.figure(figsize=(8,6))
	for i, j in z:
        	plt.plot(i[20000::1000, (1)], j[::100, ic])
	plt.xlabel('% Helicity', fontsize=24)
	plt.ylabel('tIC ' + str(ic + 1), fontsize=24)
	plt.tight_layout()
	plt.savefig(outfile + '-SvsIC' + str(ic + 1) + '.eps')
	plt.clf()

def colorWIN(data,ic,outfile):
	""" plot tIC 1 vs another IC colored by US window
	data: TICA data
	ic: give the tIC you want to plot
	outfile: name of plot
	"""
	for i in data:
		plt.plot(i[::100, 0], i[::100, ic])
	plt.xlabel('tIC 1', fontsize=18)
	plt.ylabel('tIC ' + str(ic + 1), fontsize=18)
	plt.tight_layout()
	plt.savefig(outfile + 'IC1v' + str(ic + 1) + 'colorWIN.eps')
	plt.clf()

def colorCV(z,ic,outfile):
	""" plot tIC 1 vs another IC colored by CV
	z: give zipped colvar file and tica
	ic: give the tIC you want to plot
	outfile: name of plot
	"""
	set_pub()
	fig = plt.figure(figsize=(8,6))
	for i, j in z:
		if edge == True:
			plt.scatter(i[::100, 0], i[::100, ic], c=j[20000::1000, 1], vmin=0, vmax=1, edgecolors='k', linewidths=0.25)
		else:
			plt.scatter(i[::100, 0], i[::100, ic], c=j[20000::1000, 1], vmin=0, vmax=1)
	plt.xlabel('tIC 1', fontsize=18)
	plt.ylabel('tIC ' + str(ic + 1), fontsize=24)
	cb = plt.colorbar()
	cb.set_label('% Helicity', fontsize=24)
	plt.tight_layout()
	plt.savefig(outfile + 'IC1v' + str(ic + 1) + 'colorCV.eps')
	plt.clf()

def colorTIME(z,ic,outfile):
	""" plot tIC 1 vs another IC colored by simulation time
	z: give zipped colvar file and tica
	ic: give the tIC you want to plot
	outfile: name of plot
	"""
	set_pub()
	fig = plt.figure(figsize=(8,6))
	for i, j in z:
		if edge == True:
			plt.scatter(i[::100, 0], i[::100, ic], c=(j[20000::1000, 0]/1000), vmin=20, vmax=80, edgecolors='k', linewidths=0.25)
		else:
			plt.scatter(i[::100, 0], i[::100, ic], c=(j[20000::1000, 0]/1000), vmin=20, vmax=80)
	plt.xlabel('tIC 1', fontsize=24)
	plt.ylabel('tIC ' + str(ic + 1), fontsize=24)
	cb = plt.colorbar()
	cb.set_label('time (ns)', fontsize=24)
	plt.tight_layout()
	plt.savefig(outfile + 'IC1v'+ str(ic + 1) + 'colorTIME.eps')
	plt.clf()

for i in range(0,ics):
	checklin(zipped, i, out)

for i in range(1,ics):
	plotFES(tica, i, out)
	colorWIN(tica, i, out)
	colorCV(zipped2, i, out)
	colorTIME(zipped3, i, out)
