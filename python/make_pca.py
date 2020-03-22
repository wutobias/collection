#!/usr/bin/env python

import argparse
import pytraj as pt
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from plot_hist import hist
import numpy as np
from copy import deepcopy as cp

parser = argparse.ArgumentParser(description='Process some PCA on MD trajectories.')
parser.add_argument('--traj',      '-t',  type=str, help='Trajectory file', required=True)
parser.add_argument('--parm',      '-p',  type=str, help='Parameter file', required=True)
parser.add_argument('--prefix',    '-pf', type=str, help='File name for results plot', default="PCA")
parser.add_argument('--start',     '-s',  type=int, help='First frame', default=0)
parser.add_argument('--mask',      '-m',  type=str, help='Selection Mask', default="!@H=")
parser.add_argument('--evec',      '-ev', type=int, help='Number of eigenvectors to calculate.', default=3)
parser.add_argument('--riter',     '-ri', type=int, help='Number of random trajectories to calcualte for z-score.', default=500)
parser.add_argument('--zscore',    '-z',  type=int, help='Perform z-score calculation.')
#parser.add_argument('--mask',      '-m',  type=str, help='Selection Mask', default="@CA")
#parser.add_argument('--mask',      '-m',  type=str, help='Selection Mask', default="!@11-20,H=")
#parser.add_argument('--mask',      '-m',  type=str, help='Selection Mask', default="!@1-10,H=")
parser.add_argument('--traj_proj',  '-tp', type=str, help='Additional trajectory file that is used for PC projection')
parser.add_argument('--parm_proj',  '-pp', type=str, help='Additional parameter file that is used for PC projection')
parser.add_argument('--start_proj', '-sp', type=int, help='First frame projection', default=0)
parser.add_argument('--mask_proj',  '-mp', type=str, help='Selection Mask Projection.')


def main():

	traj     = pt.load(args.traj, args.parm)

	rnd_iter = args.riter
	rnd_vecs = args.evec
	pairs    = list()

	if args.mask_proj == None:
		args.mask_proj = args.mask

	print "Mask     : ", args.mask
	print "Mask proj: ", args.mask_proj

	if rnd_vecs < 1:

		rnd_vecs = 3*traj[args.mask].xyz.shape[1] - 6

	#make pairs
	for n_i in range(rnd_vecs):
		for n_j in range(rnd_vecs):

			if n_i < n_j:

				pairs.append((n_i, n_j))

	sele     = pt.select(traj.top, args.mask)

	sele_txt = ""

	for s_i, s in enumerate(sele):

		sele_txt += "%d %s\n" %(s_i, traj.top.atomlist[s])

	o = open("%s_sele.dat" %args.prefix, "w")
	o.write(sele_txt)
	o.close()

	n_vecs          = rnd_vecs
	pca_data, eigen = pt.pca(traj[args.start:], mask=args.mask, n_vecs=n_vecs)
	eigen_val       = eigen[0]
	eigen_vec       = eigen[1]
	np.savetxt("%s_eigen_vec.dat" %args.prefix, np.c_[eigen_vec[0], eigen_vec[1], eigen_vec[2]])
	np.savetxt("%s_pcadata.dat" %args.prefix, pca_data.T)

	#h        = hist(pca_data[0], pca_data[1])
	#h.plot2d(xlab="PC1 [$\AA$]", ylab="PC2 [$\AA$]", title="PCA", name=args.out)

	# Plot PCA
	for pc_i, pc_j in pairs:

		plt.scatter(pca_data[pc_i], pca_data[pc_j], marker='o', c="r", alpha=0.5)
		plt.xlabel("PC%d [$\AA$]" %pc_i)
		plt.ylabel("PC%d [$\AA$]" %pc_j)
		plt.title("PCA PC%d vs. PC%d" %(pc_i, pc_j))
		plt.savefig("PC%d-vs-PC%s_%s.png" %(pc_i, pc_j, args.prefix))
		plt.close('all')

	# Plot atom contritbuion
	for pc_i in range(3):

		l = eigen_vec[pc_i].shape[0]
		c = np.linalg.norm(eigen_vec[pc_i].reshape((l/3, 3)), axis=1)
		a = np.arange(l/3) + 1
		plt.plot(a, c, label="PC%s"%pc_i, alpha=0.5)
		plt.legend()

	plt.xlim(0, l/3 + 1)
	plt.xlabel("Atom ID")
	plt.ylabel("Eigenvector components")
	plt.title("Eigenvectors")
	plt.savefig("Eigenvectors_%s.png" %args.prefix)
	plt.close('all')

	total_var = np.sum(eigen_val)

	plt.scatter(range(1,n_vecs+1), (np.cumsum(eigen_val)/total_var)*100, label="Cumulative Variance")
	plt.plot(range(1,n_vecs+1), (eigen_val/total_var)*100, "g--", label="Variance")
	plt.legend()
	#plt.xticks(range(1, n_vecs+1, 2))
	plt.xlabel("Eigenvector #")
	plt.ylabel("Variance explained [%]")
	plt.title("Variance explained by PC Eigenvectors")
	plt.savefig("Variance_%s.png" %args.prefix, dpi=1000)
	plt.close('all')

	if args.traj_proj != None and args.parm_proj != None:

		traj_proj       = pt.load(args.traj_proj, args.parm_proj)
		pt.rmsd(traj_proj, mask=args.mask_proj)
		#avg_proj        = pt.mean_structure(traj_proj, mask=args.mask)
		#pt.rmsd(traj_proj, mask=args.mask, ref=avg_proj)
		projection_data = pt.projection(traj_proj[args.start_proj:], args.mask_proj, eigenvalues=eigen_val,\
                                                                                     eigenvectors=eigen_vec,\
                                                                                     scalar_type='covar')
		np.savetxt("%s_pcadata_proj.dat" %args.prefix, projection_data.T)

		#h = hist(projection_data[0], projection_data[1])
		#h.plot2d(xlab="PC1 [$\AA$]", ylab="PC2 [$\AA$]", title="PCA projection", name=args.out_proj)
		for pc_i, pc_j in pairs:

			plt.scatter(pca_data[pc_i], pca_data[pc_j], marker='o', c="r", alpha=0.5)
			plt.scatter(projection_data[pc_i], projection_data[pc_j], marker='o', c="g", alpha=0.5)
			plt.xlabel("PC%d [$\AA$]" %pc_i)
			plt.ylabel("PC%d [$\AA$]" %pc_j)
			plt.title("PCA PC%d vs. PC%d with projection" %(pc_i, pc_j))
			plt.savefig("PC%d-vs-PC%s_%s_projection.png" %(pc_i, pc_j, args.prefix))
			plt.close('all')

			plt.scatter(projection_data[pc_i], projection_data[pc_j], marker='o', c="g", alpha=0.5)
			plt.xlabel("PC%d [$\AA$]" %pc_i)
			plt.ylabel("PC%d [$\AA$]" %pc_j)
			plt.title("PCA PC%d vs. PC%d only projection" %(pc_i, pc_j))
			plt.savefig("PC%d-vs-PC%d_%s_only_projection.png" %(pc_i, pc_j, args.prefix))
			plt.close('all')

		pca_data_2, eigen_2 = pt.pca(traj_proj[args.start_proj:], mask=args.mask_proj, n_vecs=n_vecs)
		eigen_val_2         = eigen_2[0]
		eigen_vec_2         = eigen_2[1]

		overlap = 0

		for pc_i in range(rnd_vecs):

			for pc_j in range(rnd_vecs):

				overlap += (np.dot(eigen_vec[pc_i], eigen_vec_2[pc_j])/(np.linalg.norm(eigen_vec[pc_i])*np.linalg.norm(eigen_vec_2[pc_j])))**2

		overlap /= rnd_vecs
		print "Vector space spanned by traj-1 overlap with traj-2 subspace (%d vecs): %6.3f" %(rnd_vecs, overlap)

		if args.zscore != None:

			overlap_rnd = np.zeros(rnd_iter)

			for r in range(rnd_iter):

				### make random traj
				t1_rnd = traj
				for f in range(t1_rnd.xyz[args.start:].shape[0]):

					idxs      = np.arange(t1_rnd.xyz[args.start+f,].shape[0])
					sele      = np.random.permutation(idxs)
					t1_rnd[f] = t1_rnd.xyz[args.start+f,][sele]

				pca_t1_rnd, eigen_t1_rnd = pt.pca(t1_rnd[args.start:], mask=args.mask, n_vecs=n_vecs)

				eigen_vec_1_rnd          = eigen_t1_rnd[1]

				for pc_i in range(n_vecs):

					for pc_j in range(n_vecs):

						overlap_rnd[r] += (np.dot(eigen_vec[pc_i], eigen_vec_1_rnd[pc_j])/(np.linalg.norm(eigen_vec[pc_i])*np.linalg.norm(eigen_vec_1_rnd[pc_j])))**2

				overlap_rnd[r] /= n_vecs

			z_score = (overlap - np.mean(overlap_rnd) ) / np.std(overlap_rnd)

			print "Z-score                                                              : %6.3f" %z_score

if __name__ == "__main__":

	args = parser.parse_args()

	main()