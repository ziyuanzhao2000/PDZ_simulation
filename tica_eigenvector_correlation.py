# written by Ziyuan Zhao, 01082022
# For a molecular system, the tICA eigenvectors associated with the k largest, distinct eigenvalues gives relative
# weights to order parameters that contributes to the k slowest-relaxing / highest autocorrelation dof. Here, the
# script reads in the tICA model calculated from two molecular systems and compute Pearson correlation coeff. between
# every pair of eigenvectors with eigenvalues above a threshold.

