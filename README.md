
# Correlation-Based Iterative Clustering for Time Course Gene Expression Data

This repository contains the MATLAB code used in the publication:

**"Correlation-based iterative clustering methods for time course data: The identification of temporal gene response modules for influenza infection in humans"**  
*Michelle Carey, Shuang Wu, Guojun Gan, Hulin Wu*  
Published in *Infectious Disease Modelling*, 2016.

## 🧬 Overview

This code implements two novel clustering algorithms designed to identify inhomogeneous gene expression modules from time course microarray data:

- **Iterative Hierarchical Clustering (IHC)**
- **Iterative Pairwise-correlation Clustering (IPC)**

These methods were applied to human gene expression data collected during an influenza infection challenge, and outperformed existing techniques such as the Markov Cluster Algorithm (MCL) and Generalized Mixture Models (GMM) in simulation and real-world settings.

## 📁 Contents

- `IHC.m` – Main function implementing the Iterative Hierarchical Clustering algorithm  
- `IPC.m` – Main function implementing the Iterative Pairwise-correlation Clustering algorithm  
- `Clustering_methods_GSE30550.m` – Script to analyze the GSE30550 influenza gene expression dataset

## 🧪 Requirements

- MATLAB R2014b or later

## 🚀 Usage

1. Preprocess or load your gene expression time course data (e.g., from GEO GSE30550).
2. Choose a correlation threshold (e.g., `alpha = 0.7`) to control cluster tightness.
3. Run the desired clustering method:
   ```matlab
   load('GSE30550.mat')

   %Subject
   sub = 1;

   %Correlation Threshold
   alpha = 0.10;

   %IHC
   [fidxcluster,rmclusters,c,mean_clusters_mat,clusters]=IHC(gexpst{sub},alpha);
   % fidxcluster: index of the genes in each cluster
   % rmclusters: clusters being removed at the prune step
   % c: convergence %of genes that are consectiviely being assigned to the same cluster
   % mean_clusters_mat: cluster centers
   % clusters

   %IPC
   [fidxcluster,rmclusters,c,mean_clusters_mat,clusters]=IPC(gexpst{sub},alpha);
   % fidxcluster: index of the genes in each cluster
   % rmclusters: clusters being removed at the prune step
   % c: convergence %of genes that are consectiviely being assigned to the same cluster
   % mean_clusters_mat: cluster centers
   % clusters
   ```

## 📖 Citation

If you use this code, please cite the original paper:

Carey M, Wu S, Gan G, Wu H. *Correlation-based iterative clustering methods for time course data: The identification of temporal gene response modules for influenza infection in humans.* Infectious Disease Modelling. 2016;1:28–39. [https://doi.org/10.1016/j.idm.2016.07.001](https://doi.org/10.1016/j.idm.2016.07.001)
