#!/usr/bin/env python3
"""
Author:   Ling Zhang & Xuewei Cao
Version:  V 0.1.0
Data:     2021.12
Describe: Tools for solving the penalized regression models, such as HuberNet, HuberLasso,
          Huber Elastic Net, Elastic Net, Lasso, Network regularization, using Accelerated
          Proximal Gradient Descent (APGD) algorithm.
Github:   https://github.com/tobefuture/APGD
Contact:  lingzhan@mtu.edu; xueweic@mtu.edu
"""

## APGD functions

## Required Modules
import pandas as pd
import numpy as np
import cvxpy as cp
from numpy import *
import networkx as nx
import time
import os
import sys
from random import sample
from sklearn.preprocessing import scale
from sklearn.linear_model import Ridge, RidgeCV
from sklearn.metrics import mean_squared_error
import warnings
warnings.filterwarnings("ignore")


## Read data by file path
## read_file(file_path)
# Read file path to get data
# input
#	file_path: the path of the file (.txt .csv) and separator by tab ('\t').
# output
#	df: data frame 
def read_file(file_path):
    file_part = pd.read_csv(file_path, sep='\t', index_col=0, chunksize=10000000)
    whole = []
    for part in file_part:
        whole.append(part)
    df = pd.concat(whole)
    print("Loaded data.")
    return df


## Simulation parts
##1.ConstructNetwork(n_genes, structure)
# Construct the network structure from either Hierarchical Network or Barabasi-Albert Network in simulation studies
# input
#	n_genes: the number of genes
#	structure: "HN": Hierarchical Network or  "BAN": Barabasi-Albert Network
# output
#	adj_all: n_genes * n_genes dimensional symmetric adjacency matrix of network structure.
def ConstructNetwork(n_genes, structure):
    ### Count n_module
    # Hierarchical Network as scenario
    if structure == "HN":
        n_each = 100
        if n_genes % n_each != 0:
            print("Error: In Hierarchical Network, number of genes must be the integer times 100!")
            sys.exit()
        else:
            n_module = n_genes / n_each
        n_module = int(n_module)
        ## Construct Hierarchical Networks
        adj = np.zeros((100, 100))
        adj[0, np.arange(1, 99, step=11)] = 1
        for ig in range(9):
            adj_group = np.zeros((11, 11))
            adj_group[0, (1, 2, 3, 4)] = 1
            adj_group[1, 5] = adj_group[2, (6, 7)] = adj_group[3, (8, 9)] = adj_group[4, 10] = 1
            start = 11 * ig + 1
            end = 11 * (ig + 1) + 1
            adj[start:end, start:end] = adj_group
        adj = adj + adj.T
        # Adjacency matrix for all network
        adj_all = np.zeros((100 * n_module, 100 * n_module))
        for im in range(n_module):
            start = 100 * im
            end = 100 * (im + 1)
            adj_all[start:end, start:end] = adj
    # Barabasi-Albert Network as scenario
    elif structure == "BAN":
        n_each = 10
        if n_genes % n_each != 0:
            print("Error: In Barabasi-Albert Network, number of genes must be the integer times 10!")
            sys.exit()
        else:
            n_module = n_genes / n_each
        n_module = int(n_module)
        # Adjacency matrix for all network for Barabasi-Albert Networks
        adj_all = np.zeros((10 * n_module, 10 * n_module))
        for im in range(n_module):
            g = nx.barabasi_albert_graph(10, 4)
            adj = nx.adjacency_matrix(g).todense()
            start = 10 * im
            end = 10 * (im + 1)
            adj_all[start:end, start:end] = adj
    # Structure error
    else:
        print(
            "Error: Please check 'structure', must be 'HN' as Hierarchical Network or 'BAN' as Barabasi-Albert Network!")
        sys.exit()
    return adj_all


##2.GraphicalModel(adj, a1=-0.7, a2=-0.1, b1=0.1, b2=0.7)
# Simulate a covariance matrix from the specific graph (Adj) based on a Gaussian graphical model.
# Input
#	adj: the adjacency matrix of network structure.
#	a1, a2, b1, b2: parameters for constructing domain [a1, a2] union [b1, b2]
#					default: a1 = -0.7, a2 = -0.1, b1 = 0.1, b2 = 0.7
# Output
#	sigma: covariance matrix of target genes according to network structure.
def GraphicalModel(adj, a1=-0.7, a2=-0.1, b1=0.1, b2=0.7):
    ## Check adjacency matrix
    if adj.shape[0] != adj.shape[1]:
        print("Error: The adjacency matrix must be squared matrix!")
        sys.exit()
    else:
        dia_vector = np.zeros(adj.shape[1])
        if not (adj == adj.T).all():
            print("Error: The adjacency matrix must be symmetric matrix!")
            sys.exit()
        elif not np.array_equal(np.diag(adj), dia_vector):
            print("Error: The diagonal elements of adjacency matrix must be all 0!")
            sys.exit()
        elif not np.array_equal(np.unique(adj), [0, 1]):
            print("Error: The elements of adjacency matrix must be 0 or 1!")
            sys.exit()
    ## Calculate Sigma
    n = adj.shape[1] ** 2
    a = a1
    b = a2
    c = b1
    d = b2
    rho = np.random.uniform(0, b - a + d - c, n)
    rho1 = (a + rho) * (rho < (b - a)) + (c + rho - (b - a)) * (rho >= (b - a))
    adj = adj * reshape(rho1, (adj.shape[1], adj.shape[1]))
    # Rescale matrix
    rsum = np.sum(abs(adj), axis=1)
    adj_rescale = adj / (5 * rsum[:, None])
    # Ensure symmetric
    sym = 0.5 * (adj_rescale + adj_rescale.T)
    np.fill_diagonal(sym, 1)
    # Calculate Sigma
    sym_inv = np.linalg.pinv(sym)
    sym_diag = mat(diag(sym_inv))
    sigma = sym_inv / sqrt(np.dot(sym_diag.T, sym_diag))
    return sigma


# 3.SimulationData(n_samples, n_genes, adj, sigma, beta0=None, beta_true=None)
# Simulate y and X from a given network structure.
# input
#	n_samples: the number of sample size
#	n_genes: the number of target genes
#	adj: the adjacency matrix of network structure. Adjacency matrix must be a n_genes * n_genes dimensional
#	symmetric matrix, the elements equal 1 indicates two genes are connected. If you consider Barabasi-Albert
#	Network or Hierarchical Network in the article, you can directly use "ConstructNetwork" function to get
#	the adjacency matrix.the adjacency matrix of network structure. directly use "ConstructNetwork" function
#	to get the adjacency matrix.
#	sigma: the covariance matrix of target genes according to network structure. You can directly use
#	"GraphicalModel" function to get the covariance matrix.
#	method: "HN": by Hierarchical Network, "BAN": by Barabasi-Albert Network or "DIY": by user designed
#	beta0: numeric value of effect size in simulation settings. 
#		   default: None; if method is "HN" or "BAN", input a numerical value.
#	beta_true: numeric matrix with the dimension of n_genes * 1 in simulation settings. 
#		   default: None; if method is "DIY", input a numerical matrix (n_genes * 1).
# output
#	y: expression levels of a transcription factor (TF)
#	X: expression levels of n_genes target genes (TGs)
#	beta: true regulated effect beta for n_genes TGs.
def SimulationData(n_samples, n_genes, adj, sigma, method, beta0=None, beta_true=None):
    ## Check method:
    # "DIY"
    if method == "DIY":
        if beta_true is not None:
            if int(np.asmatrix(beta_true).shape[0]) != n_genes:
                print("Error: In 'DIY' method, beta_true must be a numerical matrix with the dimension of n_genes * 1!")
                sys.exit()
        else:
            print(
                "Error: In 'DIY' method, please give a numerical matrix with the dimension of n_genes * 1 to beta_true!")
            sys.exit()
    # "HN" and "BAN"
    elif method == "HN" or method == "BAN":
        if beta0 is not None:
            if int(np.asmatrix(beta0).shape[0]) != 1:
                print("Error: In 'HN' and 'BAN' method, beta0 must be a numerical value (not vector or other types)!")
                sys.exit()
        else:
            print("Error: In 'HN' and 'BAN' method, please give a numerical value to beta0 (ex. 0.1)!")
            sys.exit()
    # Other typo
    else:
        print("Please check method, must be one of HN, BAN or DIY!")
        sys.exit()
    ## Check adjacency matrix
    if adj.shape[0] != adj.shape[1]:
        print("Error: The adjacency matrix must be squared matrix!")
        sys.exit()
    else:
        dia_vector = np.zeros(adj.shape[1])
        if not (adj == adj.T).all():
            print("Error: The adjacency matrix must be symmetric matrix!")
            sys.exit()
        elif not np.array_equal(np.diag(adj), dia_vector):
            print("Error: The diagonal elements of adjacency matrix must be all 0!")
            sys.exit()
        elif not np.array_equal(np.unique(adj), [0, 1]):
            print("Error: The elements of adjacency matrix must be 0 or 1!")
            sys.exit()
    ## Check covariance matrix
    if sigma.shape[0] != sigma.shape[1]:
        print("Error: The covariance matrix must be squared matrix!")
        sys.exit()
    else:
        dia_vector2 = np.ones(sigma.shape[1])
        if not (abs(sigma - sigma.T) < 1e-10).all():
            print("Error: The covariance matrix must be symmetric matrix!")
            sys.exit()
        elif not np.array_equal(np.diag(sigma), dia_vector2):
            print("Error: The diagonal elements of covariance matrix must be all 1!")
            sys.exit()

    # step 1: generate y
    y = mat(np.random.normal(0, 1, n_samples)).T
    # step 2: calculate degree
    n_degree = np.sum(adj, axis=0)
    # step 3: generate X
    beta = np.zeros((n_genes, 1))
    # Construct beta from Hierarchical Network
    if method == "HN":
        beta[0, ] = beta0
        pos1 = np.append(np.arange(1, 11 + 1), np.arange(23, 33 + 1))
        beta[pos1] = sqrt(n_degree[pos1, None]) * beta0 / 3
        pos2 = np.append(np.arange(12, 22 + 1), np.arange(34, 44 + 1))
        beta[pos2] = - sqrt(n_degree[pos2, None]) * beta0 / 3
    # Construct beta Barabasi-Albert Network
    elif method == "BAN":
        pos1 = np.append(np.arange(0, 9 + 1), np.arange(20, 30 + 1))
        beta[pos1] = sqrt(n_degree[pos1, None]) * beta0
        pos2 = np.append(np.arange(10, 20 + 1), np.arange(30, 40 + 1))
        beta[pos2] = - sqrt(n_degree[pos2, None]) * beta0
    # Construct beta beta by user matrix
    elif method == "DIY":
        beta = np.asmatrix(beta_true)
    z = np.random.multivariate_normal(np.zeros(n_genes), sigma, n_samples)
    error = np.random.multivariate_normal(np.zeros(n_genes), np.identity(n_genes), n_samples)
    X = np.dot(y, beta.T) + z + error
    return y, X, beta


## Real Data parts
##1.CalculateAdj(annotated_matrix)
# Calculate adjacency matrix from an annotation file
# input
#	annotated_matrix: n_genes * n_pathways dimensional matrix that indicates the annotation of genes
#	within pathways information.
# output
#	adj: the adjacency matrix of network structure.
def CalculateAdj(annotated_matrix):
    pw_gene_vs_pw_matrix = annotated_matrix.to_numpy()
    # Adjacency matrix
    adj = np.dot(pw_gene_vs_pw_matrix, pw_gene_vs_pw_matrix.T)
    # Set those genes have at least one same pathway, set value as 1
    adj = np.where(adj > 0, 1, 0)
    # Set diagonal value as 0
    np.fill_diagonal(adj, 0)
    pwg_name = annotated_matrix.index
    row_sum = np.sum(adj, axis=0)
    if sum(row_sum == 0) != 0:
        remove_pwg = pwg_name[where(row_sum == 0)]
        print("Please remove below unconnected %d PWGs from input files (annotation and expression files): " % (sum(row_sum == 0)) + str(
            remove_pwg))
    return adj


##2.CalculateLaplacian(adj)
# Calculate Laplacian matrix and symmetric normalized Laplacian matrix from an adjacency matrix.
# input
#	adj: the adjacency matrix of network structure.
# output
#	l: the Laplacian matrix for network structure.
#	l_norm: the symmetric normalized Laplacian matrix for network structure.
def CalculateLaplacian(adj):
    ## Check adjacency matrix
    if adj.shape[0] != adj.shape[1]:
        print("Error: The adjacency matrix must be squared matrix!")
        sys.exit()
    else:
        dia_vector = np.zeros(adj.shape[1])
        if not (adj == adj.T).all():
            print("Error: The adjacency matrix must be symmetric matrix!")
            sys.exit()
        elif not np.array_equal(np.diag(adj), dia_vector):
            print("Error: The diagonal elements of adjacency matrix must be all 0!")
            sys.exit()
        elif not np.array_equal(np.unique(adj), [0, 1]):
            print("Error: The elements of adjacency matrix must be 0 or 1!")
            sys.exit()
    # Degree matrix
    d = mat(diag(np.sum(adj, axis=1)))
    # Laplacian matrix
    l = d - adj
    l_norm = l / sqrt(diag(d))
    l_norm = l_norm.T / sqrt(diag(d))
    return l, l_norm


##3.Lambda_grid(X, y, n_lambda, alpha, loss_func)
# Simulate a grid set of lambdas for a given alpha in penalized regression.
# input
# 	X: expression levels of n_genes target genes (TGs)
# 	y: expression levels of a transcription factor (TF)
# 	n_lambda: the number of lambdas. Positive integers.
# 	alpha: the proportion of l1 norm affects (the numerical values of nonzero coefficients)
# 		   it's in range (0,1].
# 	loss_func: either "Huber" or "MSE". 
# 		  If flag = "Huber", the loss function in penalized regression model is Huber function. 
# 		  If flag = "MSE", the loss function in penalized regression model is mean squared errors.
# output
# 	lambdas: n_lambda length vector of lambdas according to the alpha you provided.
def Lambda_grid(X, y, n_lambda, alpha, loss_func):
    if type(X) is not matrix:
        X = X.to_numpy()
        X = scale(X)
    else:
        X = scale(X)
    if type(y) is not matrix:
        y = y.to_numpy()
        y = scale(y)
    else:
        y = scale(y)
    ## Check X y
    if X.shape[0] != y.shape[0]:
        print("Error: Please check the sample size of TGs and TF. They should be the same!")
        sys.exit()
    elif y.shape[1] != 1:
        print("Error: Please check the number of TF. It should be 1!")
        sys.exit()
    ## Check alpha range
    if alpha <= 0 or alpha > 1:
        print("Error: The range of alpha should be in (0,1]!")
        sys.exit()
    ratio = 0.01
    ## Check loss function
    # Use Huber function
    if loss_func == "Huber":
        m = np.subtract(*np.percentile(np.dot(X.T, y), [75, 25])) / 1.345
        h_grad = Huber_gradient(y, m)
        lambda_max = (np.abs(np.dot(X.T, h_grad)).max() / alpha)
    # Use mean squared errors function
    elif loss_func == "MSE":
        lambda_max = (np.abs(np.dot(X.T, y)).max() / alpha)
    else:
        print("Error: Please set the loss function from 'Huber' or 'MSE'!")
        sys.exit()
    if lambda_max <= np.finfo(float).resolution:
        lambdas = np.empty(n_lambda)
        lambdas.fill(np.finfo(float).resolution)
        return lambdas
    return np.logspace(np.log10(lambda_max * ratio), np.log10(lambda_max), num=n_lambda)[::-1]


##4.HuberNet_Beta(X, y, adj, lambda0, alpha0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False)
# Estimate beta_hat using HuberNet function
# input
#	X: expression levels of n_genes target genes (TGs)
#	y: expression levels of a transcription factor (TF)
#	adj: the adjacency matrix of network structure.
#	lambda0: one of parameters in HuberNet regression, which controls the number of nonzero coefficients.
#	alpha0: one of parameters in HuberNet regression, which controls the numerical values of nonzero coefficients.
#   method: The current methods must be 'APGD' or 'CVX'.
#	gamma: initial value of gamma in APGD.
#		   default: 1000
#	niter: the maximum number of APGD to solve HuberNet regression.
#		   default: 2000
#   crit_beta: converge criterion of change of beta.
#              default: 1e-4
#   crit_obj: converge criterion of change of objective function.
#             default: 1e-8
#   quiet: decide if exist the output report.
#          default: FALSE
# output
#	beta_hat: n_genes length vector of estimated regulated effect sizes, where beta_j != 0 indicates
#	j th gene is not selected in HuberNet regression.
def HuberNet_Beta(X, y, adj, lambda0, alpha0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False):
    if type(X) is not matrix:
        X = X.to_numpy()
        X = scale(X)
    else:
        X = scale(X)
    if type(y) is not matrix:
        y = y.to_numpy()
        y = scale(y)
    else:
        y = scale(y)
    ## Check X y
    if X.shape[0] != y.shape[0]:
        print("Error: Please check the sample size of TGs and TF. They should be the same!")
        sys.exit()
    elif y.shape[1] != 1:
        print("Error: Please check the number of TF. It should be 1!")
        sys.exit()
    ## Check alpha0 range
    if alpha0 <= 0 or alpha0 > 1:
        print("Error: The range of alpha should be in (0,1]!")
        sys.exit()
    ## Check lambda0
    if lambda0 <= 0:
        print("Error: Lambda should be lager than 0!")
        sys.exit()
    ## Check gamma
    if gamma <= 0:
        print("Error: Gamma should be lager than 0!")
        sys.exit()
    ## Check adjacency matrix
    if adj.shape[0] != adj.shape[1]:
        print("Error: The adjacency matrix must be squared matrix!")
        sys.exit()
    else:
        dia_vector = np.zeros(adj.shape[1])
        if not (adj == adj.T).all():
            print("Error: The adjacency matrix must be symmetric matrix!")
            sys.exit()
        elif not np.array_equal(np.diag(adj), dia_vector):
            print("Error: The diagonal elements of adjacency matrix must be all 0!")
            sys.exit()
        elif not np.array_equal(np.unique(adj), [0, 1]):
            print("Error: The elements of adjacency matrix must be 0 or 1!")
            sys.exit()
    ## Check method
    if method != "APGD" and method != "CVX":
        print("Error: The current methods must be 'APGD' or 'CVX'!")
        sys.exit()

    #  Setting:
    p = X.shape[1]
    # Laplacian matrix:
    l, l_norm = CalculateLaplacian(adj)
    # Calculate S: ridge regression:
    alpha_seq = 10 ** np.linspace(start=-2, stop=10, num=301)
    reg_cv = RidgeCV(alphas=alpha_seq, scoring='neg_mean_squared_error', normalize=True).fit(X, y)
    best_beta_hat = reg_cv.coef_
    # Sign matrix:
    s = np.zeros([p, p])
    np.fill_diagonal(s, np.sign(best_beta_hat))
    # Calculate LSL:
    sls = multiply(l_norm, diag(s))
    sls = multiply(sls.T, diag(s))
    # M: shape parameter
    m = np.subtract(*np.percentile(np.dot(X.T, y), [75, 25])) / 1.345
    if not quiet:
        print("\nStart calculating beta using Huber-Net by %s method:" % method)
    # APGD
    if method == "APGD":
        iter_num = 0
        beta_est = zeros([p, niter + 2])
        obj = np.zeros([niter + 1, 1])
        y = mat(y.T)
        for k in range(1, niter):
            # APGD python
            beta_k = beta_est[:, [k]]
            # objective function
            z = y.T - np.dot(X, beta_k)
            obj[k] = sum(Huber_Mz(z, m)) + lambda0 * alpha0 * sum(np.abs(beta_k)) + 0.5 * lambda0 * (
                    1 - alpha0) * np.dot(np.dot(beta_k.T, sls), beta_k)
            # ksi definition
            ksi = beta_k + k / (k + 3) * (beta_k - beta_est[:, [k - 1]])
            # compute the gradient
            delta = y.T - np.dot(X, ksi)
            grad_g_ksi = -np.dot(X.T, (Huber_gradient(delta, m))) + lambda0 * (1 - alpha0) * np.dot(sls, ksi)

            while True:
                theta = ksi - gamma * grad_g_ksi
                beta_prox = multiply(np.sign(theta), np.maximum(np.abs(theta) - gamma * lambda0 * alpha0, 0))
                z_prox = y.T - np.dot(X, beta_prox)
                g_beta_prox = np.sum(Huber_Mz(z_prox, m)) + 0.5 * lambda0 * (1 - alpha0) * np.dot(
                    np.dot(beta_prox.T, sls), beta_prox)
                z_ksi = y.T - np.dot(X, ksi)
                g_ksi = np.sum(Huber_Mz(z_ksi, m)) + 0.5 * lambda0 * (1 - alpha0) * np.dot(np.dot(ksi.T, sls), ksi)
                g_hat = g_ksi + np.dot(grad_g_ksi.T, (beta_prox - ksi)) + 0.5 / gamma * np.sum(
                    np.square(beta_prox - ksi))
                if g_beta_prox > g_hat:
                    gamma = 0.5 * gamma
                else:
                    break
            beta_est[:, [k + 1]] = beta_prox
            z1 = y.T - np.dot(X, beta_prox)
            obj[[k + 1]] = np.sum(Huber_Mz(z1, m)) + lambda0 * alpha0 * np.sum(abs(beta_prox)) + 0.5 * lambda0 * (1 - alpha0) * np.dot(np.dot(beta_prox.T, sls), beta_prox)
            diff_beta = beta_prox - beta_k
            temp_beta = np.sum(np.abs(diff_beta) >= crit_beta)
            temp_obj = (abs(obj[k + 1] - obj[k]) < crit_obj) * 1
            iter_num = k
            if temp_obj == 1 or temp_beta == 0:
                break
        if iter_num == niter-1:
            print("Warning: Setting number of iterations doesn't reach to the convergence. Please set larger 'niter'!")
            #sys.exit()
        beta_hat = beta_est[:, [iter_num]]
        if sum(beta_hat != 0) == 0:
            print("Warning: All estimated regression coefficients (beta_hat) are 0. Please check your size of lambda or input files!")
    # CVX
    else:
        beta = cp.Variable((p, 1))
        # try:
        #     exec("obj = cp.sum(cp.huber(y - X*beta, m)) + lambda0*alpha0*cp.pnorm(beta, 1) + 0.5*lambda0*(1-alpha0)*cp.quad_form(beta, sls)")
        # except SyntaxError:
        obj = cp.sum(cp.huber(y - X @ beta, m)) + lambda0 * alpha0 * cp.pnorm(beta, 1) + 0.5 * lambda0 * (1 - alpha0) * cp.quad_form(beta, sls)
        prob = cp.Problem(cp.Minimize(obj))
        results = prob.solve()
        beta_hat = beta.value
    return beta_hat


##5.HuberNet_SP(X, y, adj, alphas, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True)
# Estimate selection probability using HuberNet function solving by APGD
# input
#	X: expression levels of n_genes target genes (TGs)
#	y: expression levels of a transcription factor (TF)
#   adj: the adjacency matrix of network structure.
#	alphas: the grid sets of alpha (in [0,1]) used to calculate selection probabilities of genes.
#	n_lambda: the number of lambdas
#	B: the number of half-sample resampling used to calculate selection probabilities of genes.
#	   default: 500
#	gamma: initial value of gamma in APGD.
#		   default: 1000
#	niter: the maximum number of APGD to solve HuberNet regression.
#		   default: 2000
#   crit_beta: converge criterion of change of beta.
#              default: 1e-4
#   crit_obj: converge criterion of change of objective function.
#             default: 1e-8
#   timer: decide if exist the output report.
#          default: True
# output
#	sp_hubernet: n_genes length vector of selection probability.
def HuberNet_SP(X, y, adj, alphas, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True):
    X_ori = X
    y_ori = y
    if type(X) is not matrix:
        X = X.to_numpy()
        X = scale(X)
    else:
        X = scale(X)
    if type(y) is not matrix:
        y = y.to_numpy()
        y = scale(y)
    else:
        y = scale(y)
    ## Check X y
    if X.shape[0] != y.shape[0]:
        print("Error: Please check the sample size of TGs and TF. They should be the same!")
        sys.exit()
    elif y.shape[1] != 1:
        print("Error: Please check the number of TF. It should be 1!")
        sys.exit()
    ## Check adjacency matrix
    if adj.shape[0] != adj.shape[1]:
        print("Error: The adjacency matrix must be squared matrix!")
        sys.exit()
    else:
        dia_vector = np.zeros(adj.shape[1])
        if not (adj == adj.T).all():
            print("Error: The adjacency matrix must be symmetric matrix!")
            sys.exit()
        elif not np.array_equal(np.diag(adj), dia_vector):
            print("Error: The diagonal elements of adjacency matrix must be all 0!")
            sys.exit()
        elif not np.array_equal(np.unique(adj), [0, 1]):
            print("Error: The elements of adjacency matrix must be 0 or 1!")
            sys.exit()
    ## Check alphas
    for alpha in alphas:
        if alpha <= 0 or alpha > 1:
            print("Error: The range of alpha should be in (0,1]!")
            sys.exit()

    #  Setting:
    p = X.shape[1]
    n = X.shape[0]

    ## Start Huber-Net by APGD method
    print("\nStart calculating selection probability using Huber-Net by APGD method:")
    # Set the cases
    num_alpha_set = len(alphas)
    sp_la = empty([p, n_lambda * num_alpha_set])  # Selection probability matrix 2539 * cases
    temp = 0
    # See the cases results of SP values and beta_hat values.
    start_apgd = time.time()
    for j in range(len(alphas)):
        alpha = alphas[j]
        lamda_set = Lambda_grid(X_ori, y_ori, n_lambda, alpha, "Huber")
        for i in range(len(lamda_set)):
            lamda = lamda_set[i]
            temp += 1
            print("\nCase " + str(temp) + ": lambda " + str(lamda) + ", alpha " + str(alpha))
            beta_hat_la = empty([p, B])
            prop_la = empty([p, B])
            # Resampling B times:
            start_case = time.time()
            for b in range(B):
                sub_pos = sample(list(range(n)), int(n / 2))
                X_sub = matrix(X[sub_pos, :])
                y_sub = matrix(y[sub_pos, :])
                # APGD python
                beta_hat_la[:, [b]] = HuberNet_Beta(X_sub, y_sub, adj, lamda, alpha, method="APGD", gamma=gamma, niter=niter, crit_beta=crit_beta, crit_obj=crit_obj, quiet=True)
            end_case = time.time()
            if timer:
                count_time(start_case, end_case)
            prop_la = np.where(beta_hat_la != 0, 1, 0)
            sp_la[:, [temp - 1]] = mat(prop_la.mean(axis=1)).T
    sp_hubernet = np.max(sp_la, axis=1)

    # save_file() function
    # savetxt(output_path+str(index)+'_SP_hubernet_'+TF_name+'.csv', SP_hubernet, delimiter=',')

    print("\nDone with Huber-Net by APGD method!")
    end_apgd = time.time()
    if timer:
        count_time(start_apgd, end_apgd)
    return sp_hubernet


##6.HuberLasso_Beta(X, y, lambda0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False)
# Estimate beta_hat using Huber Lasso function
# input
#	X: expression levels of n_genes target genes (TGs)
#	y: expression levels of a transcription factor (TF)
#	lambda0: one of parameters in Huber Lasso regression, which controls the number of nonzero coefficients.
#   method: the current methods must be 'APGD' or 'CVX'.
#	gamma: initial value of gamma in APGD.
#		   default: 1000
#	niter: the maximum number of APGD to solve Huber Lasso regression.
#		   default: 2000
#   crit_beta: converge criterion of change of beta.
#              default: 1e-4
#   crit_obj: converge criterion of change of objective function.
#             default: 1e-8
#   quiet: decide if exist the output report.
#          default: FALSE
# output
#	beta_hat: n_genes length vector of estimated regulated effect sizes, where beta_j != 0 indicates
#	j th gene is not selected in HuberLasso regression.
def HuberLasso_Beta(X, y, lambda0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False):
    if type(X) is not matrix:
        X = X.to_numpy()
        X = scale(X)
    else:
        X = scale(X)
    if type(y) is not matrix:
        y = y.to_numpy()
        y = scale(y)
    else:
        y = scale(y)
    ## Check X y
    if X.shape[0] != y.shape[0]:
        print("Error: Please check the sample size of TGs and TF. They should be the same!")
        sys.exit()
    elif y.shape[1] != 1:
        print("Error: Please check the number of TF. It should be 1!")
        sys.exit()
    ## Check lambda0
    if lambda0 <= 0:
        print("Error: Lambda should be lager than 0!")
        sys.exit()
    ## Check gamma
    if gamma <= 0:
        print("Error: Gamma should be lager than 0!")
        sys.exit()
    ## Check method
    if method != "APGD" and method != "CVX":
        print("Error: The current methods must be 'APGD' or 'CVX'!")
        sys.exit()

    #  Setting:m
    p = X.shape[1]
    # m: shape parameter
    m = np.subtract(*np.percentile(np.dot(X.T, y), [75, 25])) / 1.345
    if not quiet:
        print("\nStart calculating beta using Huber-Lasso by %s method:" % method)
    # APGD
    if method == "APGD":
        iter_num = 0
        beta_est = zeros([p, niter + 2])
        obj = np.zeros([niter + 1, 1])
        y = mat(y.T)
        for k in range(1, niter):
            # APGD python
            beta_k = beta_est[:, [k]]
            # objective function
            z = y.T - np.dot(X, beta_k)
            obj[k] = sum(Huber_Mz(z, m)) + lambda0 * sum(np.abs(beta_k))
            # ksi definition
            ksi = beta_k + k / (k + 3) * (beta_k - beta_est[:, [k - 1]])
            # compute the gradient
            delta = y.T - np.dot(X, ksi)
            grad_g_ksi = -np.dot(X.T, (Huber_gradient(delta, m)))

            while True:
                theta = ksi - gamma * grad_g_ksi
                beta_prox = multiply(np.sign(theta), np.maximum(np.abs(theta) - gamma * lambda0, 0))
                z_prox = y.T - np.dot(X, beta_prox)
                g_beta_prox = np.sum(Huber_Mz(z_prox, m))
                z_ksi = y.T - np.dot(X, ksi)
                g_ksi = np.sum(Huber_Mz(z_ksi, m))
                g_hat = g_ksi + np.dot(grad_g_ksi.T, (beta_prox - ksi)) + 0.5 / gamma * np.sum(
                    np.square(beta_prox - ksi))
                if g_beta_prox > g_hat:
                    gamma = 0.5 * gamma
                else:
                    break
            beta_est[:, [k + 1]] = beta_prox
            z1 = y.T - np.dot(X, beta_prox)
            obj[[k + 1]] = np.sum(Huber_Mz(z1, m)) + lambda0 * np.sum(abs(beta_prox))
            diff_beta = beta_prox - beta_k
            temp_beta = np.sum(np.abs(diff_beta) >= crit_beta)
            temp_obj = (abs(obj[k + 1] - obj[k]) < crit_obj) * 1
            iter_num = k
            if temp_obj == 1 or temp_beta == 0:
                break
        if iter_num == niter - 1:
            print("Warning: Setting number of iterations doesn't reach to the convergence. Please set larger 'niter'!")
            # sys.exit()
        beta_hat = beta_est[:, [iter_num]]
        if sum(beta_hat != 0) == 0:
            print("Warning: All estimated regression coefficients (beta_hat) are 0. Please check your size of lambda or input files!")
    # CVX
    else:
        beta = cp.Variable((p, 1))
        # try:
        #     exec("cp.sum(cp.huber(y - X@beta, M)) + lambda0*cp.pnorm(beta, 1)")
        # except SyntaxError:
        obj = cp.sum(cp.huber(y - X @ beta, m)) + lambda0 * cp.pnorm(beta, 1)
        prob = cp.Problem(cp.Minimize(obj))
        results = prob.solve()
        beta_hat = beta.value

    return beta_hat


##7.HuberLasso_SP(X, y, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True)
# Estimate selection probability using HuberLasso function solving by APGD
# input
#	X: expression levels of n_genes target genes (TGs)
#	y: expression levels of a transcription factor (TF)
#	alphas: the grid sets of alpha (in [0,1]) used to calculate selection probabilities of genes.
#	n_lambda: the number of lambdas
#	B: the number of half-sample resampling used to calculate selection probabilities of genes.
#	   default: 500
#	gamma: initial value of gamma in APGD.
#		   default: 1000
#	niter: the maximum number of APGD to solve Huber Lasso regression.
#		   default: 2000
#   crit_beta: converge criterion of change of beta.
#              default: 1e-4
#   crit_obj: converge criterion of change of objective function.
#             default: 1e-8
#   timer: decide if exist the output report.
#          default: True
# output
#	sp_huberlasso: n_genes length vector of selection probability.
def HuberLasso_SP(X, y, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True):
    X_ori = X
    y_ori = y
    if type(X) is not matrix:
        X = X.to_numpy()
        X = scale(X)
    else:
        X = scale(X)
    if type(y) is not matrix:
        y = y.to_numpy()
        y = scale(y)
    else:
        y = scale(y)
    ## Check X y
    if X.shape[0] != y.shape[0]:
        print("Error: Please check the sample size of TGs and TF. They should be the same!")
        sys.exit()
    elif y.shape[1] != 1:
        print("Error: Please check the number of TF. It should be 1!")
        sys.exit()

    #  Setting:
    p = X.shape[1]
    n = X.shape[0]

    ## Start Huber-Net by APGD method
    print("\nStart calculating selection probability using Huber-Lasso by APGD method:")
    # TF_name = y.columns
    # print("TF gene: "+TF_name)
    # Set the cases
    lamda_set = Lambda_grid(X_ori, y_ori, n_lambda, 1, "Huber")
    print(lamda_set)
    sp_la = empty([p, len(lamda_set)])
    temp = 0
    # See the cases results of SP values and beta_hat values.
    start_apgd = time.time()
    for i in range(len(lamda_set)):
        lamda = lamda_set[i]
        temp += 1
        print("\nCase " + str(temp) + ": lambda " + str(lamda))
        beta_hat_la = empty([p, B])
        prop_la = empty([p, B])
        # Resampling B times:
        start_case = time.time()
        for b in range(B):
            sub_pos = sample(list(range(n)), int(n / 2))
            X_sub = matrix(X[sub_pos, :])
            y_sub = matrix(y[sub_pos, :])
            # APGD python
            beta_hat_la[:, [b]] = HuberLasso_Beta(X_sub, y_sub, lamda, method="APGD", gamma=gamma, niter=niter, crit_beta=crit_beta, crit_obj=crit_obj, quiet=True)
        end_case = time.time()
        if timer:
            count_time(start_case, end_case)
        prop_la = np.where(beta_hat_la != 0, 1, 0)
        sp_la[:, [temp - 1]] = mat(prop_la.mean(axis=1)).T
    sp_huberlasso = np.max(sp_la, axis=1)

    # save_file() function
    # savetxt(output_path+str(index)+'_SP_huberlasso_'+TF_name+'.csv', SP_huberlasso, delimiter=',')

    print("\nDone with Huber-Lasso by APGD method!")
    end_apgd = time.time()
    if timer:
        count_time(start_apgd, end_apgd)
    return sp_huberlasso


##8.HuberEnet_Beta(X, y, lambda0, alpha0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False)
# Estimate beta_hat using Huber Elastic Net function
# input
#	X: expression levels of n_genes target genes (TGs)
#	y: expression levels of a transcription factor (TF)
#	lambda0: one of parameters in Huber Elastic Net regression, which controls the number of nonzero coefficients.
#	alpha0: one of parameters in Huber Elastic Net regression, which controls the numerical values of nonzero coefficients.
#   method: the current methods must be 'APGD' or 'CVX'.
#	gamma: initial value of gamma in APGD.
#		   default: 1000
#	niter: the maximum number of APGD to solve Huber Elastic Net regression.
#		   default: 2000
#   crit_beta: converge criterion of change of beta.
#              default: 1e-4
#   crit_obj: converge criterion of change of objective function.
#             default: 1e-8
#   quiet: decide if exist the output report.
#          default: FALSE
# output
#	beta_hat: n_genes length vector of estimated regulated effect sizes, where beta_j != 0 indicates
#	j th gene is not selected in Huber Elastic Net regression.
def HuberEnet_Beta(X, y, lambda0, alpha0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False):
    if type(X) is not matrix:
        X = X.to_numpy()
        X = scale(X)
    else:
        X = scale(X)
    if type(y) is not matrix:
        y = y.to_numpy()
        y = scale(y)
    else:
        y = scale(y)
    ## Check X y
    if X.shape[0] != y.shape[0]:
        print("Error: Please check the sample size of TGs and TF. They should be the same!")
        sys.exit()
    elif y.shape[1] != 1:
        print("Error: Please check the number of TF. It should be 1!")
        sys.exit()
    ## Check alpha0 range
    if alpha0 <= 0 or alpha0 > 1:
        print("Error: The range of alpha should be in (0,1]!")
        sys.exit()
    ## Check lambda0
    if lambda0 <= 0:
        print("Error: Lambda should be lager than 0!")
        sys.exit()
    ## Check gamma
    if gamma <= 0:
        print("Error: Gamma should be lager than 0!")
        sys.exit()
    ## Check method
    if method != "APGD" and method != "CVX":
        print("Error: The current methods must be 'APGD' or 'CVX'!")
        sys.exit()

    #  Setting:
    p = X.shape[1]
    # M: shape parameter
    m = np.subtract(*np.percentile(np.dot(X.T, y), [75, 25])) / 1.345
    if not quiet:
        print("\nStart calculating beta using Huber-Enet by %s method:" % method)
    # APGD
    if method == "APGD":
        iter_num = 0
        beta_est = zeros([p, niter + 2])
        obj = np.zeros([niter + 1, 1])
        y = mat(y.T)
        for k in range(1, niter):
            # APGD python
            beta_k = beta_est[:, [k]]
            # objective function
            z = y.T - np.dot(X, beta_k)
            obj[k] = sum(Huber_Mz(z, m)) + lambda0 * alpha0 * sum(np.abs(beta_k)) + 0.5 * lambda0 * (
                    1 - alpha0) * np.dot(beta_k.T, beta_k)
            # ksi definition
            ksi = beta_k + k / (k + 3) * (beta_k - beta_est[:, [k - 1]])
            # compute the gradient
            delta = y.T - np.dot(X, ksi)
            grad_g_ksi = -np.dot(X.T, (Huber_gradient(delta, m))) + lambda0 * (1 - alpha0) * ksi

            while True:
                theta = ksi - gamma * grad_g_ksi
                beta_prox = multiply(np.sign(theta), np.maximum(np.abs(theta) - gamma * lambda0 * alpha0, 0))
                z_prox = y.T - np.dot(X, beta_prox)
                g_beta_prox = np.sum(Huber_Mz(z_prox, m)) + 0.5 * lambda0 * (1 - alpha0) * np.dot(beta_prox.T,
                                                                                                  beta_prox)
                z_ksi = y.T - np.dot(X, ksi)
                g_ksi = np.sum(Huber_Mz(z_ksi, m)) + 0.5 * lambda0 * (1 - alpha0) * np.dot(ksi.T, ksi)
                g_hat = g_ksi + np.dot(grad_g_ksi.T, (beta_prox - ksi)) + 0.5 / gamma * np.sum(
                    np.square(beta_prox - ksi))
                if g_beta_prox > g_hat:
                    gamma = 0.5 * gamma
                else:
                    break
            beta_est[:, [k + 1]] = beta_prox
            z1 = y.T - np.dot(X, beta_prox)
            obj[[k + 1]] = np.sum(Huber_Mz(z1, m)) + lambda0 * alpha0 * np.sum(abs(beta_prox)) + 0.5 * lambda0 * (
                    1 - alpha0) * np.dot(beta_prox.T, beta_prox)
            diff_beta = beta_prox - beta_k
            temp_beta = np.sum(np.abs(diff_beta) >= crit_beta)
            temp_obj = (abs(obj[k + 1] - obj[k]) < crit_obj) * 1
            iter_num = k
            if temp_obj == 1 or temp_beta == 0:
                break
        if iter_num == niter - 1:
            print("Warning: Setting number of iterations doesn't reach to the convergence. Please set larger 'niter'!")
            # sys.exit()
        beta_hat = beta_est[:, [iter_num]]
        if sum(beta_hat != 0) == 0:
            print("Warning: All estimated regression coefficients (beta_hat) are 0. Please check your size of lambda or input files!")
    # CVX
    else:
        beta = cp.Variable((p, 1))
        # try:
        #     exec(
        #         "obj = cp.sum(cp.huber(y - X@beta, m)) + lambda0*alpha0*cp.pnorm(beta, 1) + 0.5*lambda0*(1-alpha0)*cp.pnorm(beta, p=2)**2")
        # except SyntaxError:
        obj = cp.sum(cp.huber(y - X @ beta, m)) + lambda0 * alpha0 * cp.pnorm(beta, 1) + 0.5 * lambda0 * (1 - alpha0) * cp.pnorm(beta, p=2) ** 2
        prob = cp.Problem(cp.Minimize(obj))
        results = prob.solve()
        beta_hat = beta.value
    return beta_hat


##9.HuberEnet_SP(X, y, alphas, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True)
# Estimate selection probability using HuberENET function solving by APGD
# input
#	X: expression levels of n_genes target genes (TGs)
#	y: expression levels of a transcription factor (TF)
#	alphas: the grid sets of alpha (in [0,1]) used to calculate selection probabilities of genes.
#	n_lambda: the number of lambdas
#	B: the number of half-sample resampling used to calculate selection probabilities of genes.
#	   default: 500
#	gamma: initial value of gamma in APGD.
#		   default: 1000
#	niter: the maximum number of APGD to solve Huber Elastic Net regression.
#		   default: 2000
#   crit_beta: converge criterion of change of beta.
#              default: 1e-4
#   crit_obj: converge criterion of change of objective function.
#             default: 1e-8
#   timer: decide if exist the output report.
#          default: True
# output
#	sp_huberenet: n_genes length vector of selection probability.
def HuberEnet_SP(X, y, alphas, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True):
    X_ori = X
    y_ori = y
    if type(X) is not matrix:
        X = X.to_numpy()
        X = scale(X)
    else:
        X = scale(X)
    if type(y) is not matrix:
        y = y.to_numpy()
        y = scale(y)
    else:
        y = scale(y)
    ## Check X y
    if X.shape[0] != y.shape[0]:
        print("Error: Please check the sample size of TGs and TF. They should be the same!")
        sys.exit()
    elif y.shape[1] != 1:
        print("Error: Please check the number of TF. It should be 1!")
        sys.exit()
    ## Check alphas
    for alpha in alphas:
        if alpha <= 0 or alpha > 1:
            print("Error: The range of alpha should be in (0,1]!")
            sys.exit()

    #  Setting:
    p = X.shape[1]
    n = X.shape[0]

    ## Start Huber-Net by APGD method
    print("\nStart calculating selection probability using Huber-Enet by APGD method:")
    # TF_name = y.columns
    # print("TF gene: "+TF_name)
    # Set the cases
    num_alpha_set = len(alphas)
    sp_la = empty([p, n_lambda * num_alpha_set])  # Selection probability matrix
    temp = 0
    # See the case results of SP values and beta_hat values.
    start_apgd = time.time()
    for j in range(len(alphas)):
        alpha = alphas[j]
        lamda_set = Lambda_grid(X_ori, y_ori, n_lambda, alpha, "Huber")
        for i in range(len(lamda_set)):
            lamda = lamda_set[i]
            temp += 1
            print("\nCase " + str(temp) + ": lambda " + str(lamda) + ", alpha " + str(alpha))
            beta_hat_la = empty([p, B])
            prop_la = empty([p, B])
            # Resampling B times:
            start_case = time.time()
            for b in range(B):
                sub_pos = sample(list(range(n)), int(n / 2))
                X_sub = matrix(X[sub_pos, :])
                y_sub = matrix(y[sub_pos, :])
                # APGD python
                beta_hat_la[:, [b]] = HuberEnet_Beta(X_sub, y_sub, lamda, alpha, method="APGD", gamma=gamma, niter=niter, crit_beta=crit_beta, crit_obj=crit_obj, quiet=True)
            end_case = time.time()
            if timer:
                count_time(start_case, end_case)
            prop_la = np.where(beta_hat_la != 0, 1, 0)
            sp_la[:, [temp - 1]] = mat(prop_la.mean(axis=1)).T
    sp_huberenet = np.max(sp_la, axis=1)

    # save_file() function
    # savetxt(output_path+str(index)+'_SP_hubernet_'+TF_name+'.csv', SP_huberenet, delimiter=',')

    print("\nDone with Huber-Enet by APGD method!")
    end_apgd = time.time()
    if timer:
        count_time(start_apgd, end_apgd)
    return sp_huberenet


##10.Enet_Beta(X, y, lambda0, alpha0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False)
# Estimate beta_hat using Elastic Net function
# input
#	X: expression levels of n_genes target genes (TGs)
#	y: expression levels of a transcription factor (TF)
#	lambda0: one of parameters in Elastic Net regression, which controls the number of nonzero coefficients.
#	alpha0: one of parameters in Elastic Net regression, which controls the numerical values of nonzero coefficients.
#   method: the current methods must be 'APGD' or 'CVX'.
#	gamma: initial value of gamma in APGD.
#		   default: 1000
#	niter: the maximum number of APGD to solve Elastic Net regression.
#		   default: 2000
#   crit_beta: converge criterion of change of beta.
#              default: 1e-4
#   crit_obj: converge criterion of change of objective function.
#             default: 1e-8
#   quiet: decide if exist the output report.
#          default: FALSE
# output
#	beta_hat: n_genes length vector of estimated regulated effect sizes, where beta_j != 0 indicates
#	j th gene is not selected in Elastic Net regression.
def Enet_Beta(X, y, lambda0, alpha0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False):
    if type(X) is not matrix:
        X = X.to_numpy()
        X = scale(X)
    else:
        X = scale(X)
    if type(y) is not matrix:
        y = y.to_numpy()
        y = scale(y)
    else:
        y = scale(y)
    ## Check X y
    if X.shape[0] != y.shape[0]:
        print("Error: Please check the sample size of TGs and TF. They should be the same!")
        sys.exit()
    elif y.shape[1] != 1:
        print("Error: Please check the number of TF. It should be 1!")
        sys.exit()
    ## Check alpha0 range
    if alpha0 <= 0 or alpha0 > 1:
        print("Error: The range of alpha should be in (0,1]!")
        sys.exit()
    ## Check lambda0
    if lambda0 <= 0:
        print("Error: Lambda should be lager than 0!")
        sys.exit()
    ## Check gamma
    if gamma <= 0:
        print("Error: Gamma should be lager than 0!")
        sys.exit()
    ## Check method
    if method != "APGD" and method != "CVX":
        print("Error: The current methods must be 'APGD' or 'CVX'!")
        sys.exit()

    #  Setting:
    p = X.shape[1]
    if not quiet:
        print("\nStart calculating beta using Enet by %s method:" % method)
    # APGD
    if method == "APGD":
        iter_num = 0
        beta_est = zeros([p, niter + 2])
        obj = np.zeros([niter + 1, 1])
        y = mat(y.T)
        for k in range(1, niter):
            # APGD python
            beta_k = beta_est[:, [k]]
            # objective function
            z = y.T - np.dot(X, beta_k)
            obj[k] = 0.5 * sum(np.square(z)) + lambda0 * alpha0 * sum(np.abs(beta_k)) + 0.5 * lambda0 * (
                    1 - alpha0) * np.dot(beta_k.T, beta_k)
            # ksi definition
            ksi = beta_k + k / (k + 3) * (beta_k - beta_est[:, [k - 1]])
            # compute the gradient
            delta = y.T - np.dot(X, ksi)
            grad_g_ksi = -np.dot(X.T, delta) + lambda0 * (1 - alpha0) * ksi

            while True:
                theta = ksi - gamma * grad_g_ksi
                beta_prox = multiply(np.sign(theta), np.maximum(np.abs(theta) - gamma * lambda0 * alpha0, 0))
                z_prox = y.T - np.dot(X, beta_prox)
                g_beta_prox = np.sum(np.square(z_prox)) + 0.5 * lambda0 * (1 - alpha0) * np.dot(beta_prox.T, beta_prox)
                z_ksi = y.T - np.dot(X, ksi)
                g_ksi = np.sum(np.square(z_ksi)) + 0.5 * lambda0 * (1 - alpha0) * np.dot(ksi.T, ksi)
                g_hat = g_ksi + np.dot(grad_g_ksi.T, (beta_prox - ksi)) + 0.5 / gamma * np.sum(
                    np.square(beta_prox - ksi))
                if g_beta_prox > g_hat:
                    gamma = 0.5 * gamma
                else:
                    break
            beta_est[:, [k + 1]] = beta_prox
            z1 = y.T - np.dot(X, beta_prox)
            obj[[k + 1]] = 0.5 * np.sum(np.square(z1)) + lambda0 * alpha0 * np.sum(abs(beta_prox)) + 0.5 * lambda0 * (
                    1 - alpha0) * np.dot(beta_prox.T, beta_prox)
            diff_beta = beta_prox - beta_k
            temp_beta = np.sum(np.abs(diff_beta) >= crit_beta)
            temp_obj = (abs(obj[k + 1] - obj[k]) < crit_obj) * 1
            iter_num = k
            if temp_obj == 1 or temp_beta == 0:
                break
        if iter_num == niter - 1:
            print("Warning: Setting number of iterations doesn't reach to the convergence. Please set larger 'niter'!")
            # sys.exit()
        beta_hat = beta_est[:, [iter_num]]
        if sum(beta_hat != 0) == 0:
            print("Warning: All estimated regression coefficients (beta_hat) are 0. Please check your size of lambda or input files!")
    # CVX
    else:
        beta = cp.Variable((p, 1))
        # try:
        #     exec(
        #         "obj = cp.pnorm(y - X@beta, p=2)**2 + lambda0*alpha0*cp.pnorm(beta, 1) + 0.5*lambda0*(1-alpha0)*cp.pnorm(beta, p=2)**2")
        # except SyntaxError:
        obj = cp.pnorm(y - X @ beta, p=2) ** 2 + lambda0 * alpha0 * cp.pnorm(beta, 1) + 0.5 * lambda0 * (1 - alpha0) * cp.pnorm(beta, p=2) ** 2
        prob = cp.Problem(cp.Minimize(obj))
        results = prob.solve()
        beta_hat = beta.value
    return beta_hat


##9.Enet_SP(X, y, alphas, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True)
# Estimate selection probability using ENET function solving by APGD
# input
#	X: expression levels of n_genes target genes (TGs)
#	y: expression levels of a transcription factor (TF)
#	alphas: the grid sets of alpha (in [0,1]) used to calculate selection probabilities of genes.
#	n_lambda: the number of lambdas
#	B: the number of half-sample resampling used to calculate selection probabilities of genes.
#	   default: 500
#	gamma: initial value of gamma in APGD.
#		   default: 1000
#	niter: the maximum number of APGD to solve Elastic Net regression.
#		   default: 2000
#   crit_beta: converge criterion of change of beta.
#              default: 1e-4
#   crit_obj: converge criterion of change of objective function.
#             default: 1e-8
#   timer: decide if exist the output report.
#          default: True
# output
#	sp_enet: n_genes length vector of selection probability.
def Enet_SP(X, y, alphas, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True):
    X_ori = X
    y_ori = y
    if type(X) is not matrix:
        X = X.to_numpy()
        X = scale(X)
    else:
        X = scale(X)
    if type(y) is not matrix:
        y = y.to_numpy()
        y = scale(y)
    else:
        y = scale(y)
    ## Check X y
    if X.shape[0] != y.shape[0]:
        print("Error: Please check the sample size of TGs and TF. They should be the same!")
        sys.exit()
    elif y.shape[1] != 1:
        print("Error: Please check the number of TF. It should be 1!")
        sys.exit()
    ## Check alphas
    for alpha in alphas:
        if alpha <= 0 or alpha > 1:
            print("Error: The range of alpha should be in (0,1]!")
            sys.exit()

    #  Setting:
    p = X.shape[1]
    n = X.shape[0]

    ## Start Huber-Net by APGD method
    print("\nStart calculating selection probability using Enet by APGD method:")
    # TF_name = y.columns
    # print("TF gene: "+TF_name)
    # Set the cases
    num_alpha_set = len(alphas)
    sp_la = empty([p, n_lambda * num_alpha_set])  # Selection probability matrix
    temp = 0
    # See the cases results of SP values and beta_hat values.
    start_apgd = time.time()
    for j in range(len(alphas)):
        alpha = alphas[j]
        lamda_set = Lambda_grid(X_ori, y_ori, n_lambda, alpha, "MSE")
        for i in range(len(lamda_set)):
            lamda = lamda_set[i]
            temp += 1
            print("\nCase " + str(temp) + ": lambda " + str(lamda) + ", alpha " + str(alpha))
            beta_hat_la = empty([p, B])
            prop_la = empty([p, B])
            # Resampling B times:
            start_case = time.time()
            for b in range(B):
                sub_pos = sample(list(range(n)), int(n / 2))
                X_sub = matrix(X[sub_pos, :])
                y_sub = matrix(y[sub_pos, :])
                # APGD python
                beta_hat_la[:, [b]] = Enet_Beta(X_sub, y_sub, lamda, alpha, method="APGD", gamma=gamma, niter=niter, crit_beta=crit_beta, crit_obj=crit_obj, quiet=True)
            end_case = time.time()
            if timer:
                count_time(start_case, end_case)
            prop_la = np.where(beta_hat_la != 0, 1, 0)
            sp_la[:, [temp - 1]] = mat(prop_la.mean(axis=1)).T
    sp_enet = np.max(sp_la, axis=1)

    # save_file() function
    # savetxt(output_path+str(index)+'_SP_enet_'+TF_name+'.csv', SP_enet, delimiter=',')

    print("\nDone with Enet by APGD method!")
    end_apgd = time.time()
    if timer:
        count_time(start_apgd, end_apgd)
    return sp_enet


##10.Lasso_Beta(X, y, lambda0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False)
# Estimate beta_hat using Lasso function
# input
#	X: expression levels of n_genes target genes (TGs)
#	y: expression levels of a transcription factor (TF)
#	lambda0: one of parameters in HuberNet regression, which controls the number of nonzero coefficients.
#	gamma: initial value of gamma in APGD.
#		   default: 1000
#	niter: the maximum number of APGD to solve Lasso regression.
#		   default: 2000
#   crit_beta: converge criterion of change of beta.
#              default: 1e-4
#   crit_obj: converge criterion of change of objective function.
#             default: 1e-8
#   quiet: decide if exist the output report.
#          default: FALSE
# output
#	beta_hat: n_genes length vector of estimated regulated effect sizes, where beta_j != 0 indicates
#	j th gene is not selected in Lasso regression.
def Lasso_Beta(X, y, lambda0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False):
    if type(X) is not matrix:
        X = X.to_numpy()
        X = scale(X)
    else:
        X = scale(X)
    if type(y) is not matrix:
        y = y.to_numpy()
        y = scale(y)
    else:
        y = scale(y)
    ## Check X y
    if X.shape[0] != y.shape[0]:
        print("Error: Please check the sample size of TGs and TF. They should be the same!")
        sys.exit()
    elif y.shape[1] != 1:
        print("Error: Please check the number of TF. It should be 1!")
        sys.exit()
    ## Check lambda0
    if lambda0 <= 0:
        print("Error: Lambda should be lager than 0!")
        sys.exit()
    ## Check gamma
    if gamma <= 0:
        print("Error: Gamma should be lager than 0!")
        sys.exit()
    ## Check method
    if method != "APGD" and method != "CVX":
        print("Error: The current methods must be 'APGD' or 'CVX'!")
        sys.exit()

    #  Setting:
    p = X.shape[1]
    if not quiet:
        print("\nStart calculating beta using Lasso by %s method:" % method)
    # APGD
    if method == "APGD":
        iter_num = 0
        beta_est = zeros([p, niter + 2])
        obj = np.zeros([niter + 1, 1])
        y = mat(y.T)
        for k in range(1, niter):
            # APGD python
            beta_k = beta_est[:, [k]]
            # objective function
            z = y.T - np.dot(X, beta_k)
            obj[k] = 0.5 * sum(np.square(z)) + lambda0 * sum(np.abs(beta_k))
            # ksi definition
            ksi = beta_k + k / (k + 3) * (beta_k - beta_est[:, [k - 1]])
            # compute the gradient
            delta = y.T - np.dot(X, ksi)
            grad_g_ksi = -np.dot(X.T, delta)

            while True:
                theta = ksi - gamma * grad_g_ksi
                beta_prox = multiply(np.sign(theta), np.maximum(np.abs(theta) - gamma * lambda0, 0))
                z_prox = y.T - np.dot(X, beta_prox)
                g_beta_prox = np.sum(np.square(z_prox))
                z_ksi = y.T - np.dot(X, ksi)
                g_ksi = np.sum(np.square(z_ksi))
                g_hat = g_ksi + np.dot(grad_g_ksi.T, (beta_prox - ksi)) + 0.5 / gamma * np.sum(
                    np.square(beta_prox - ksi))
                if g_beta_prox > g_hat:
                    gamma = 0.5 * gamma
                else:
                    break
            beta_est[:, [k + 1]] = beta_prox
            z1 = y.T - np.dot(X, beta_prox)
            obj[[k + 1]] = 0.5 * np.sum(np.square(z1)) + lambda0 * np.sum(abs(beta_prox))
            diff_beta = beta_prox - beta_k
            temp_beta = np.sum(np.abs(diff_beta) >= crit_beta)
            temp_obj = (abs(obj[k + 1] - obj[k]) < crit_obj) * 1
            iter_num = k
            if temp_obj == 1 or temp_beta == 0:
                break
        if iter_num == niter - 1:
            print("Warning: Setting number of iterations doesn't reach to the convergence. Please set larger 'niter'!")
            # sys.exit()
        beta_hat = beta_est[:, [iter_num]]
        if sum(beta_hat != 0) == 0:
            print("Warning: All estimated regression coefficients (beta_hat) are 0. Please check your size of lambda or input files!")
    # CVX
    else:
        beta = cp.Variable((p, 1))
        # try:
        #     exec("obj = cp.pnorm(y - X@beta, p=2)**2 + lambda0*cp.pnorm(beta, 1)")
        # except SyntaxError:
        obj = cp.pnorm(y - X @ beta, p=2) ** 2 + lambda0 * cp.pnorm(beta, 1)
        prob = cp.Problem(cp.Minimize(obj))
        results = prob.solve()
        beta_hat = beta.value
    return beta_hat


##11.Lasso_SP(X, y, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True)
# Estimate selection probability using Lasso function solving by APGD
# input
#	X: expression levels of n_genes target genes (TGs)
#	y: expression levels of a transcription factor (TF)
#	alphas: the grid sets of alpha (in [0,1]) used to calculate selection probabilities of genes.
#	n_lambda: the number of lambdas
#	B: the number of half-sample resampling used to calculate selection probabilities of genes.
#	   default: 500
#	gamma: initial value of gamma in APGD.
#		   default: 1000
#	niter: the maximum number of APGD to solve HuberNet regression.
#		   default: 2000
#   crit_beta: converge criterion of change of beta.
#              default: 1e-4
#   crit_obj: converge criterion of change of objective function.
#             default: 1e-8
#   timer: decide if exist the output report.
#          default: True
# output
#	sp_lasso: n_genes length vector of selection probability.
def Lasso_SP(X, y, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True):
    X_ori = X
    y_ori = y
    if type(X) is not matrix:
        X = X.to_numpy()
        X = scale(X)
    else:
        X = scale(X)
    if type(y) is not matrix:
        y = y.to_numpy()
        y = scale(y)
    else:
        y = scale(y)
    ## Check X y
    if X.shape[0] != y.shape[0]:
        print("Error: Please check the sample size of TGs and TF. They should be the same!")
        sys.exit()
    elif y.shape[1] != 1:
        print("Error: Please check the number of TF. It should be 1!")
        sys.exit()

    #  Setting:
    p = X.shape[1]
    n = X.shape[0]

    ## Start Huber-Net by APGD method
    print("\nStart calculating selection probability using Lasso by APGD method:")
    # TF_name = y.columns
    # print("TF gene: "+TF_name)
    # Set the cases
    lamda_set = Lambda_grid(X_ori, y_ori, n_lambda, 1, "MSE")
    print(lamda_set)
    sp_la = empty([p, len(lamda_set)])
    temp = 0
    # See the cases results of SP values and beta_hat values.
    start_apgd = time.time()
    for i in range(len(lamda_set)):
        lamda = lamda_set[i]
        temp += 1
        print("\nCase " + str(temp) + ": lambda " + str(lamda))
        beta_hat_la = empty([p, B])
        prop_la = empty([p, B])
        # Resampling B times:
        start_case = time.time()
        for b in range(B):
            sub_pos = sample(list(range(n)), int(n / 2))
            X_sub = matrix(X[sub_pos, :])
            y_sub = matrix(y[sub_pos, :])
            # APGD python
            beta_hat_la[:, [b]] = Lasso_Beta(X_sub, y_sub, lamda, method="APGD", gamma=gamma, niter=niter, crit_beta=crit_beta, crit_obj=crit_obj, quiet=True)
        end_case = time.time()
        if timer:
            count_time(start_case, end_case)
        prop_la = np.where(beta_hat_la != 0, 1, 0)
        sp_la[:, [temp - 1]] = mat(prop_la.mean(axis=1)).T
    sp_lasso = np.max(sp_la, axis=1)

    # save_file() function
    # savetxt(output_path+str(index)+'_SP_lasso_'+TF_name+'.csv', SP_lasso, delimiter=',')

    print("\nDone with Lasso by APGD method!")
    end_apgd = time.time()
    if timer:
        count_time(start_apgd, end_apgd)
    return sp_lasso


##12.Net_Beta(X, y, adj, lambda0, alpha0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False)
# Estimate beta_hat using Net function
# input
#	X: expression levels of n_genes target genes (TGs)
#	y: expression levels of a transcription factor (TF)
#	adj: the adjacency matrix of network structure.
#	lambda0: one of parameters in Net regression, which controls the number of nonzero coefficients.
#	alpha0: one of parameters in Net regression, which controls the numerical values of nonzero coefficients.
#   method: the current methods must be 'APGD' or 'CVX'.
#	gamma: initial value of gamma in APGD.
#		   default: 1000
#	niter: the maximum number of APGD to solve HuberNet regression.
#		   default: 2000
#   crit_beta: converge criterion of change of beta.
#              default: 1e-4
#   crit_obj: converge criterion of change of objective function.
#             default: 1e-8
#   quiet: decide if exist the output report.
#          default: FALSE
# output
#	beta_hat: n_genes length vector of estimated regulated effect sizes, where beta_j != 0 indicates
#	j th gene is not selected in Net regression.
def Net_Beta(X, y, adj, lambda0, alpha0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False):
    if type(X) is not matrix:
        X = X.to_numpy()
        X = scale(X)
    else:
        X = scale(X)
    if type(y) is not matrix:
        y = y.to_numpy()
        y = scale(y)
    else:
        y = scale(y)
    ## Check X y
    if X.shape[0] != y.shape[0]:
        print("Error: Please check the sample size of TGs and TF. They should be the same!")
        sys.exit()
    elif y.shape[1] != 1:
        print("Error: Please check the number of TF. It should be 1!")
        sys.exit()
    ## Check alpha0 range
    if alpha0 <= 0 or alpha0 > 1:
        print("Error: The range of alpha should be in (0,1]!")
        sys.exit()
    ## Check lambda0
    if lambda0 <= 0:
        print("Error: Lambda should be lager than 0!")
        sys.exit()
    ## Check gamma
    if gamma <= 0:
        print("Error: Gamma should be lager than 0!")
        sys.exit()
    ## Check adjacency matrix
    if adj.shape[0] != adj.shape[1]:
        print("Error: The adjacency matrix must be squared matrix!")
        sys.exit()
    else:
        dia_vector = np.zeros(adj.shape[1])
        if not (adj == adj.T).all():
            print("Error: The adjacency matrix must be symmetric matrix!")
            sys.exit()
        elif not np.array_equal(np.diag(adj), dia_vector):
            print("Error: The diagonal elements of adjacency matrix must be all 0!")
            sys.exit()
        elif not np.array_equal(np.unique(adj), [0, 1]):
            print("Error: The elements of adjacency matrix must be 0 or 1!")
            sys.exit()
    ## Check method
    if method != "APGD" and method != "CVX":
        print("Error: The current methods must be 'APGD' or 'CVX'!")
        sys.exit()

    #  Setting:
    p = X.shape[1]
    # Laplacian matrix:
    l, l_norm = CalculateLaplacian(adj)
    # Calculate S: ridge regression:
    alpha_seq = 10 ** np.linspace(start=-2, stop=10, num=301)
    reg_cv = RidgeCV(alphas=alpha_seq, scoring='neg_mean_squared_error', normalize=True).fit(X, y)
    best_beta_hat = reg_cv.coef_
    # Sign matrix:
    s = np.zeros([p, p])
    np.fill_diagonal(s, np.sign(best_beta_hat))
    # Calculate LSL:
    sls = multiply(l_norm, diag(s))
    sls = multiply(sls.T, diag(s))

    if not quiet:
        print("\nStart calculating beta using Net by %s method:" % method)
    # APGD
    if method == "APGD":
        iter_num = 0
        beta_est = zeros([p, niter + 2])
        obj = np.zeros([niter + 1, 1])
        y = mat(y.T)
        for k in range(1, niter):
            # APGD python
            beta_k = beta_est[:, [k]]
            # objective function
            z = y.T - np.dot(X, beta_k)
            obj[k] = 0.5 * sum(np.square(z)) + lambda0 * alpha0 * sum(np.abs(beta_k)) + 0.5 * lambda0 * (
                    1 - alpha0) * np.dot(np.dot(beta_k.T, sls), beta_k)
            # ksi definition
            ksi = beta_k + k / (k + 3) * (beta_k - beta_est[:, [k - 1]])
            # compute the gradient
            delta = y.T - np.dot(X, ksi)
            grad_g_ksi = -np.dot(X.T, delta) + lambda0 * (1 - alpha0) * np.dot(sls, ksi)

            while True:
                theta = ksi - gamma * grad_g_ksi
                beta_prox = multiply(np.sign(theta), np.maximum(np.abs(theta) - gamma * lambda0 * alpha0, 0))
                z_prox = y.T - np.dot(X, beta_prox)
                g_beta_prox = np.sum(np.square(z_prox)) + 0.5 * lambda0 * (1 - alpha0) * np.dot(
                    np.dot(beta_prox.T, sls), beta_prox)
                z_ksi = y.T - np.dot(X, ksi)
                g_ksi = np.sum(np.square(z_ksi)) + 0.5 * lambda0 * (1 - alpha0) * np.dot(np.dot(ksi.T, sls), ksi)
                g_hat = g_ksi + np.dot(grad_g_ksi.T, (beta_prox - ksi)) + 0.5 / gamma * np.sum(
                    np.square(beta_prox - ksi))
                if g_beta_prox > g_hat:
                    gamma = 0.5 * gamma
                else:
                    break
            beta_est[:, [k + 1]] = beta_prox
            z1 = y.T - np.dot(X, beta_prox)
            obj[[k + 1]] = 0.5 * np.sum(np.square(z1)) + lambda0 * alpha0 * np.sum(abs(beta_prox)) + 0.5 * lambda0 * (
                    1 - alpha0) * np.dot(np.dot(beta_prox.T, sls), beta_prox)
            diff_beta = beta_prox - beta_k
            temp_beta = np.sum(np.abs(diff_beta) >= crit_beta)
            temp_obj = (abs(obj[k + 1] - obj[k]) < crit_obj) * 1
            iter_num = k
            if temp_obj == 1 or temp_beta == 0:
                break
        if iter_num == niter - 1:
            print("Warning: Setting number of iterations doesn't reach to the convergence. Please set larger 'niter'!")
            # sys.exit()
        beta_hat = beta_est[:, [iter_num]]
        if sum(beta_hat != 0) == 0:
            print("Warning: All estimated regression coefficients (beta_hat) are 0. Please check your size of lambda or input files!")
    # CVX
    else:
        beta = cp.Variable((p, 1))
        # try:
        #     exec(
        #         "obj = cp.pnorm(y - X@beta, p=2)**2 + lambda0*alpha0*cp.pnorm(beta, 1) + 0.5*lambda0*(1-alpha0)*cp.quad_form(beta, SLS)")
        # except SyntaxError:
        obj = cp.pnorm(y - X @ beta, p=2) ** 2 + lambda0 * alpha0 * cp.pnorm(beta, 1) + 0.5 * lambda0 * (1 - alpha0) * cp.quad_form(beta, sls)
        prob = cp.Problem(cp.Minimize(obj))
        results = prob.solve()
        beta_hat = beta.value
    return beta_hat


##13.Net_SP(X, y, adj, alphas, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True)
# Estimate selection probability using Net function solving by APGD
# input
#	X: expression levels of n_genes target genes (TGs)
#	y: expression levels of a transcription factor (TF)
#	alphas: the grid sets of alpha (in [0,1]) used to calculate selection probabilities of genes.
#	n_lambda: the number of lambdas
#	B: the number of half-sample resampling used to calculate selection probabilities of genes.
#	   default: 500
#	gamma: initial value of gamma in APGD.
#		   default: 1000
#	niter: the maximum number of APGD to solve Net regression.
#		   default: 2000
#   crit_beta: converge criterion of change of beta.
#              default: 1e-4
#   crit_obj: converge criterion of change of objective function.
#             default: 1e-8
#   timer: decide if exist the output report.
#          default: True
# output
#	sp_net: n_genes length vector of selection probability.
def Net_SP(X, y, adj, alphas, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True):
    X_ori = X
    y_ori = y
    if type(X) is not matrix:
        X = X.to_numpy()
        X = scale(X)
    else:
        X = scale(X)
    if type(y) is not matrix:
        y = y.to_numpy()
        y = scale(y)
    else:
        y = scale(y)
    ## Check X y
    if X.shape[0] != y.shape[0]:
        print("Error: Please check the sample size of TGs and TF. They should be the same!")
        sys.exit()
    elif y.shape[1] != 1:
        print("Error: Please check the number of TF. It should be 1!")
        sys.exit()
    ## Check adjacency matrix
    if adj.shape[0] != adj.shape[1]:
        print("Error: The adjacency matrix must be squared matrix!")
        sys.exit()
    else:
        dia_vector = np.zeros(adj.shape[1])
        if not (adj == adj.T).all():
            print("Error: The adjacency matrix must be symmetric matrix!")
            sys.exit()
        elif not np.array_equal(np.diag(adj), dia_vector):
            print("Error: The diagonal elements of adjacency matrix must be all 0!")
            sys.exit()
        elif not np.array_equal(np.unique(adj), [0, 1]):
            print("Error: The elements of adjacency matrix must be 0 or 1!")
            sys.exit()
    ## Check alphas
    for alpha in alphas:
        if alpha <= 0 or alpha > 1:
            print("Error: The range of alpha should be in (0,1]!")
            sys.exit()

    #  Setting:
    p = X.shape[1]
    n = X.shape[0]

    ## Start Huber-Net by APGD method
    print("\nStart calculating selection probability using Net by APGD method:")
    # TF_name = y.columns
    # print("TF gene: "+TF_name)
    # Set the cases
    num_alpha_set = len(alphas)
    sp_la = empty([p, n_lambda * num_alpha_set])  # Selection probability matrix 2539 * cases
    temp = 0
    # See the cases results of SP values and beta_hat values.
    start_apgd = time.time()
    for j in range(len(alphas)):
        alpha = alphas[j]
        lamda_set = Lambda_grid(X_ori, y_ori, n_lambda, alpha, "MSE")
        for i in range(len(lamda_set)):
            lamda = lamda_set[i]
            temp += 1
            print("\nCase " + str(temp) + ": lambda " + str(lamda) + ", alpha " + str(alpha))
            beta_hat_la = empty([p, B])
            prop_la = empty([p, B])
            # Resampling B times:
            start_case = time.time()
            for b in range(B):
                sub_pos = sample(list(range(n)), int(n / 2))
                X_sub = matrix(X[sub_pos, :])
                y_sub = matrix(y[sub_pos, :])
                # APGD python
                beta_hat_la[:, [b]] = HuberNet_Beta(X_sub, y_sub, adj, lamda, alpha, method="APGD", gamma=gamma, niter=niter, crit_beta=crit_beta, crit_obj=crit_obj, quiet=True)
            end_case = time.time()
            if timer:
                count_time(start_case, end_case)
            prop_la = np.where(beta_hat_la != 0, 1, 0)
            sp_la[:, [temp - 1]] = mat(prop_la.mean(axis=1)).T
    sp_net = np.max(sp_la, axis=1)

    # save_file() function
    # savetxt(output_path+str(index)+'_SP_net_'+TF_name+'.csv', SP_net, delimiter=',')

    print("\nDone with Net by APGD method!")
    end_apgd = time.time()
    if timer:
        count_time(start_apgd, end_apgd)
    return sp_net


## Support function:
## Huber_gradient(delta, m)
# Calculate the gradient of Huber function for an input value delta.
# input
#   delta: Input value delta
#   m: Shape parameter, which is defaulted to be one-tenth of the interquartile range (IRQ).
# output
#   value: The gradient of Huber function for an input value z
def Huber_gradient(delta, m):
    value = np.where(np.abs(delta) <= m, 2 * delta, 2 * m * np.sign(delta))
    return value


## Huber_Mz(z, m)
# Calculate the Huber function for an input value z.
# input
#   z: Input value z
#   m: Shape parameter, which is defaulted to be one-tenth of the interquartile range (IRQ).
# output
#   value: The Huber function for an input value z
def Huber_Mz(z, m):
    value = np.where(np.abs(z) <= m, np.square(z), 2 * m * np.abs(z) - m ** 2)  # 369*1
    return value


# Timer
# Calculate running time
def count_time(start, end):
    cost_h = (end - start) / 3600
    cost_min = (end - start) / 60
    cost_s = end - start
    print("The calculate time: " + str(cost_h) + "h, " + str(cost_min) + "min, " + str(cost_s) + "s.")
    return
