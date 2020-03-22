#!/usr/bin/env python

from __future__ import print_function

import argparse

from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

import pytraj as pt

parser = argparse.ArgumentParser(description='K-means clustering on cartesian coordinates of MD trajectories.')
parser.add_argument('--traj',      '-t',  type=str, help='Trajectory file', required=True)
parser.add_argument('--parm',      '-p',  type=str, help='Parameter file', required=True)
parser.add_argument('--prefix',    '-pf', type=str, help='Prefix for file names.', default="clusteranalysis")
parser.add_argument('--mask',      '-m',  type=str, help='Selection Mask', default="@CA")
parser.add_argument('--stride',    '-st', type=int, help='Stride frames', default=10)
parser.add_argument('--start',     '-s',  type=int, help='Start at this frame.', default=0)
parser.add_argument('--pca',       '-pc', type=str, help='Perform PCA prior to Clustering.', default="no", choices=["yes", "no"])

def main():

    X     = pt.load(args.traj, args.parm, stride=args.stride)

    if args.pca == "no":

        X     = X[args.mask].xyz[args.start:]
        shape = X.shape
        X     = X.reshape((shape[0], shape[1]*3))

    else:

        n_vecs          = 3*X[args.mask].xyz[args.start:].shape[1] - 6
        pca_data, eigen = pt.pca(X[args.start:], n_vecs=n_vecs, mask=args.mask)
        eigen_val       = eigen[0]
        eigen_vec       = eigen[1]
        np.savetxt("%s_eigen_vec.dat" %args.prefix, np.c_[eigen_vec[0], eigen_vec[1], eigen_vec[2]])

        pairs  = list()

        #make pairs
        for n_i in range(3):
            for n_j in range(3):
        
                if n_i < n_j:
        
                    pairs.append((n_i, n_j))

        # Plot PCA
        for pc_i, pc_j in pairs:
    
            plt.scatter(pca_data[pc_i], pca_data[pc_j], marker='o', c="r", alpha=0.5)
            plt.xlabel("PC%d [$\AA$]" %pc_i)
            plt.ylabel("PC%d [$\AA$]" %pc_j)
            plt.title("PCA PC%d vs. PC%d" %(pc_i, pc_j))
            plt.savefig("PC%d-vs-PC%s_%s.png" %(pc_i, pc_j, args.prefix), dpi=1000)
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
        plt.savefig("Eigenvectors_%s.png" %args.prefix, dpi=1000)
        plt.close('all')
    
        total_var = np.sum(eigen_val)
    
        plt.scatter(range(1,n_vecs+1), (np.cumsum(eigen_val)/total_var)*100, label="Cumulative Variance")
        plt.plot(range(1,n_vecs+1), (eigen_val/total_var)*100, "g--", label="Eigenvector Variance")
        plt.legend()
        #plt.xticks(range(1, n_vecs+1, 2))
        plt.xlabel("Eigenvector #")
        plt.ylabel("Fractional of Variance explained [%]")
        plt.title("Explained total variance explained by PCA")
        plt.savefig("Variance_%s.png" %args.prefix, dpi=1000)
        plt.close('all')

        X = pca_data
    
    range_n_clusters = range(2,20)
    
    for n_clusters in range_n_clusters:
        # Create a subplot with 1 row and 2 columns
        if args.pca == "yes":
            fig, (ax1, ax2) = plt.subplots(1, 2)
            fig.set_size_inches(18, 7)
        else:
            fig, (ax1, ax2) = plt.subplots(1, 1)
    
        # The 1st subplot is the silhouette plot
        # The silhouette coefficient can range from -1, 1 but in this example all
        # lie within [-0.1, 1]
        ax1.set_xlim([-0.1, 1])
        # The (n_clusters+1)*10 is for inserting blank space between silhouette
        # plots of individual clusters, to demarcate them clearly.
        ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])
    
        # Initialize the clusterer with n_clusters value and a random generator
        # seed of rand for reproducibility.
        rand = np.random.randint(99999)
        print("Random seed is %d." %rand)
        clusterer = KMeans(n_clusters=n_clusters, random_state=rand)
        cluster_labels = clusterer.fit_predict(X)
    
        # The silhouette_score gives the average value for all the samples.
        # This gives a perspective into the density and separation of the formed
        # clusters
        silhouette_avg = silhouette_score(X, cluster_labels)
        print("For n_clusters =", n_clusters,
              "The average silhouette_score is :", silhouette_avg)
    
        # Compute the silhouette scores for each sample
        sample_silhouette_values = silhouette_samples(X, cluster_labels)
    
        y_lower = 10
        for i in range(n_clusters):
            # Aggregate the silhouette scores for samples belonging to
            # cluster i, and sort them
            ith_cluster_silhouette_values = \
                sample_silhouette_values[cluster_labels == i]
    
            ith_cluster_silhouette_values.sort()
    
            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i
    
            color = cm.spectral(float(i) / n_clusters)
            ax1.fill_betweenx(np.arange(y_lower, y_upper),
                              0, ith_cluster_silhouette_values,
                              facecolor=color, edgecolor=color, alpha=0.7)
    
            # Label the silhouette plots with their cluster numbers at the middle
            ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
    
            # Compute the new y_lower for next plot
            y_lower = y_upper + 10  # 10 for the 0 samples
    
        ax1.set_title("The silhouette plot for the various clusters.")
        ax1.set_xlabel("The silhouette coefficient values")
        ax1.set_ylabel("Cluster label")
    
        # The vertical line for average silhouette score of all the values
        ax1.axvline(x=silhouette_avg, color="red", linestyle="--")
    
        ax1.set_yticks([])  # Clear the yaxis labels / ticks
        ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
    
        if args.pca == "yes":

            # 2nd Plot showing the actual clusters formed
            colors = cm.spectral(cluster_labels.astype(float) / n_clusters)
            ax2.scatter(X[:, 0], X[:, 1], marker='.', s=30, lw=0, alpha=0.7,
                        c=colors)
        
            # Labeling the clusters
            centers = clusterer.cluster_centers_
            # Draw white circles at cluster centers
            ax2.scatter(centers[:, 0], centers[:, 1],
                        marker='o', c="white", alpha=1, s=200)
        
            for i, c in enumerate(centers):
                ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1, s=50)
    
            ax2.set_title("The visualization of the clustered data.")
            ax2.set_xlabel("Feature space for the 1st feature")
            ax2.set_ylabel("Feature space for the 2nd feature")
    
        plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
                      "with n_clusters = %d" % n_clusters),
                     fontsize=14, fontweight='bold')
    
        plt.savefig("%s_silhouette_n=%d.png" %(args.prefix, n_clusters), dpi=1000)
        plt.close('all')

if __name__ == "__main__":

    args = parser.parse_args()

    main()