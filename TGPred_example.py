import TGPred as tp
from numpy import *
import numpy as np
from sklearn.preprocessing import scale

adj1 = tp.ConstructNetwork(200, "BAN")
#print(adj1)
adj2 = tp.ConstructNetwork(200, "HN")
# print(adj2.shape)
#
sigma1 = tp.GraphicalModel(adj1)
#print(sigma1)
sigma2 = tp.GraphicalModel(adj2)

## test wrong ex.
# w_adj1 = np.array([[0, 1, 1],
#                    [1, 1, 1],
#                    [1, 1, 0]])
# w_sigma = tp.GraphicalModel(w_adj1)

# w_sigma = np.array([[1, 1, 1],
#                     [1, 1, 0],
#                     [1, 0, 1]])
y, X, beta = tp.SimulationData(300, 200, adj2, sigma2, "HN", beta0=10, beta_true=None)
# print(type(y))
#print(y.shape, X.shape, beta.shape)
# print(y, X, beta)
# lam = tp.Lambda_grid(X, y, 10, 0.1, "MSE")
# print(lam)

PWs_matrix_path = "./Sample_data/Matrix_131PWGs_vs_5PWs.txt"
Annoted_df = tp.read_file(PWs_matrix_path)
Adj = tp.CalculateAdj(Annoted_df)
#print(Adj, Adj.shape)
l, l_norm = tp.CalculateLaplacian(Adj)
print(l, l_norm)

#
TF_exp_path = "./Sample_data/1TF_genes_expression.txt"
PW_exp_path = "./Sample_data/131PW_genes_expression.txt"
PWs_matrix_path = "./Sample_data/Matrix_131PWGs_vs_5PWs.txt"
y = tp.read_file(TF_exp_path)
X = tp.read_file(PW_exp_path)
df = tp.read_file(PWs_matrix_path)
Adj = tp.CalculateAdj(df)
X = tp.read_file(PW_exp_path)
y = tp.read_file(TF_exp_path)

#
# TF_exp = tp.read_file(TF_exp_path)	#739*1525
# TF_genes = TF_exp.columns
# name = str(TF_genes[0])
# y = TF_exp[[name]]
#
# PW_exp = tp.read_file(PW_exp_path) #739*2539

# Hubernet
beta_hat = tp.HuberNet_Beta(X, y, Adj, 20, 0.1, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False, scales=False)
print(beta_hat.shape)
print(beta_hat)
sp = tp.HuberNet_SP(X, y, Adj, [0.1, 0.9], 2, B=4, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True)
print(sp)
print(where(sp>0))

# HuberLasso
beta_hat = tp.HuberLasso_Beta(X, y, 10, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False)
print(beta_hat.shape)
print(beta_hat)
sp = tp.HuberLasso_SP(X, y, 5, B=4, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True)
print(sp)
print(where(sp>0))

# HuberENET
beta_hat = tp.HuberENET_Beta(X, y, 10, 0.2, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False, scales=False)
print(beta_hat.shape)
print(beta_hat)
sp = tp.HuberENET_SP(X, y, [0.1,0.9], 2, B=4)
print(sp)
print(where(sp>0))

# MSENet
beta_hat = tp.MSENet_Beta(X, y, Adj, 10, 0.2, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False, scales=False)
print(beta_hat.shape)
print(beta_hat)
sp = tp.MSENet_SP(X, y, Adj, [0.1, 0.2], 2, B=4)
print(sp)
print(where(sp>0))

# MSELasso
beta_hat = tp.MSELasso_Beta(X, y, 20, method="CVX", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False, scales=False)
print(beta_hat.shape)
print(beta_hat)
sp = tp.MSELasso_SP(X, y, 2, B=4)
print(sp)
print(where(sp==1))

# MSEENET
beta_hat = tp.MSEENET_Beta(X, y, 20, 0.2, method="CVX", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False, scales=False)
print(beta_hat.shape)
print(beta_hat)
sp = tp.MSEENET_SP(X, y, [0.1, 0.2], 2, B=4, niter=3000)
print(sp)
print(where(sp==1))