<!-- README.md is generated from README.Rmd. Please edit that file -->

# APGD v.0.1.0

<!-- badges: start -->

Python version of the Accelerated Proximal Gradient Descent (APGD) algorithm is to solve the penalized regression models, including 

- **HuberNet**: Huber loss function along with Network-based penalty function;
- **HuberLasso**: Huber loss function along with Lasso penalty function;
- **HuberENET**: Huber loss function along with Elastic Net penalty function;
- **ENET**: Mean square error loss function along with Elastic Net penalty function;
- **Lasso**: Mean square error loss function along with Lasso penalty function;
- **Net**: Mean square error loss function along with Network-based penalty function.

## Functions

### 1. read_file(file_path)
**Read file path to get data.**

- Input:   
	- file_path: the path of the file (.txt .csv) and separator by tab ('\t'). 
  
- Output:   	
	- df: data frame.    

&emsp; &emsp;
   	
### 2. ConstructNetwork(n_genes, structure)
**Construct the network structure from either Hierarchical Network or Barabasi-Albert Network in simulation studies.**   

- Input:
	- n_genes: the number of genes.   
	- structure: "HN": Hierarchical Network or  "BAN": Barabasi-Albert Network.   

- Output:   
	- adj_all: n_genes * n_genes dimensional symmetric adjacency matrix of network structure.	  
	
&emsp; &emsp;

### 3. GraphicalModel(adj, a1=-0.7, a2=-0.1, b1=0.1, b2=0.7)
**Simulate a covariance matrix from the specific graph (Adj) based on a Gaussian graphical model.**  

- Input:
	- adj: the adjacency matrix of network structure.
	- a1, a2, b1, b2: parameters for constructing domain [a1, a2] union [b1, b2].  
		- default: a1 = -0.7, a2 = -0.1, b1 = 0.1, b2 = 0.7   
		
- Output: 
	- sigma: covariance matrix of target genes according to network structure.  

&emsp; &emsp;

### 4. SimulationData(n_samples, n_genes, adj, sigma, beta0=None, beta_true=None)
**Simulate y and X from a given network structure.**  

- Input:
	- n_samples: the number of sample size.
	- n_genes: the number of target genes.
	- adj: the adjacency matrix of network structure. Adjacency matrix must be a n_genes * n_genes dimensional.
	- symmetric matrix, the elements equal 1 indicates two genes are connected. If you consider Barabasi-Albert.
	- Network or Hierarchical Network in the article, you can directly use "ConstructNetwork" function to get the adjacency matrix.the adjacency matrix of network structure. directly use "ConstructNetwork" function to get the adjacency matrix.
	- sigma: the covariance matrix of target genes according to network structure. You can directly use "GraphicalModel" function to get the covariance matrix.
	- method: "HN": by Hierarchical Network, "BAN": by Barabasi-Albert Network or "DIY": by user designed.
	- beta0: numeric value of effect size in simulation settings.  
		- default: None; if method is "HN" or "BAN", input a numerical value.
	- beta_true: numeric matrix with the dimension of n_genes * 1 in simulation settings.  
		- default: None; if method is "DIY", input a numerical matrix (n_genes * 1).  
		
- Output:
	- y: expression levels of a transcription factor (TF).
	- X: expression levels of n_genes target genes (TGs).
	- beta: true regulated effect beta for n_genes TGs.  

&emsp; &emsp;

### 5. CalculateAdj(annotated_matrix)
***Calculate adjacency matrix from an annotation file.***  

- Input:
	- annotated_matrix: n_genes * n_pathways dimensional matrix that indicates the annotation of genes within pathways information.  
	
- Output:
	- adj: the adjacency matrix of network structure.  

&emsp; &emsp;

### 6. CalculateLaplacian(adj)
***Calculate Laplacian matrix and symmetric normalized Laplacian matrix from an adjacency matrix.***  

- Input:
	adj: the adjacency matrix of network structure.  
	
- Output:
	- l: the Laplacian matrix for network structure.
	- l_norm: the symmetric normalized Laplacian matrix for network structure.  

&emsp; &emsp;

### 7. Lambda_grid(X, y, n_lambda, alpha, loss_func)
***Simulate a grid set of lambdas for a given alpha in penalized regression.*** 

- Input:
	- X: expression levels of n_genes target genes (TGs).
	- y: expression levels of a transcription factor (TF).
	- n_lambda: the number of lambdas. Positive integers.
	- alpha: the proportion of l1 norm affects (the numerical values of nonzero coefficients), it's in range (0,1].
	- loss_func: either "Huber" or "MSE". 
		- If flag = "Huber", the loss function in penalized regression model is Huber function. 
		- If flag = "MSE", the loss function in penalized regression model is mean squared errors.  
		
- Output:
	- lambdas: n_lambda length vector of lambdas according to the alpha you provided.  

&emsp; &emsp;
	
### 8. HuberNet_Beta(X, y, adj, lambda0, alpha0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False)
***Estimate beta_hat using HuberNet function.***  

- Input:
	- X: expression levels of n_genes target genes (TGs).
	- y: expression levels of a transcription factor (TF).
	- adj: the adjacency matrix of network structure.
	- lambda0: one of parameters in HuberNet regression, which controls the number of nonzero coefficients.
	- alpha0: one of parameters in HuberNet regression, which controls the numerical values of nonzero coefficients.
 	- method: The current methods must be 'APGD' or 'CVX'.
	- gamma: initial value of gamma in APGD.
		- default: 1000
	- niter: the maximum number of APGD to solve HuberNet regression.
		- default: 2000
 	- crit_beta: converge criterion of change of beta.
 	 	- default: 1e-4
 	- crit_obj: converge criterion of change of objective function.
		- default: 1e-8
	- quiet: decide if exist the output report.
		- default: FALSE  
		
- Output:
	- beta_hat: n_genes length vector of estimated regulated effect sizes, where beta_j != 0 indicates j th gene is not selected in HuberNet regression.  

&emsp; &emsp;

### 9. HuberNet_SP(X, y, adj, alphas, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True)
***Estimate selection probability using HuberNet function solving by APGD.***

- Input:
	- X: expression levels of n_genes target genes (TGs).
	- y: expression levels of a transcription factor (TF).
 	- adj: the adjacency matrix of network structure.
	- alphas: the grid sets of alpha (in [0,1]) used to calculate selection probabilities of genes.
	- n_lambda: the number of lambdas
	- B: the number of half-sample resampling used to calculate selection probabilities of genes.
		- default: 500
	- gamma: initial value of gamma in APGD.
		- default: 1000
	- niter: the maximum number of APGD to solve HuberNet regression.
		- default: 2000
	- crit_beta: converge criterion of change of beta.
		- default: 1e-4
	- crit_obj: converge criterion of change of objective function.
		- default: 1e-8
	- timer: decide if exist the output report.
		- default: True  
		
- Output:
	- sp_hubernet: n_genes length vector of selection probability.  

&emsp; &emsp;

### 10. HuberLasso_Beta(X, y, lambda0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False)
***Estimate beta_hat using Huber Lasso function.***

- Input:
	- X: expression levels of n_genes target genes (TGs).
	- y: expression levels of a transcription factor (TF).
	- lambda0: one of parameters in Huber Lasso regression, which controls the number of nonzero coefficients.
	- method: the current methods must be 'APGD' or 'CVX'.
	- gamma: initial value of gamma in APGD.
		- default: 1000
	- niter: the maximum number of APGD to solve Huber Lasso regression.
		- default: 2000
 	- crit_beta: converge criterion of change of beta.
 	 	- default: 1e-4
 	- crit_obj: converge criterion of change of objective function.
		- default: 1e-8
	- quiet: decide if exist the output report.
		- default: FALSE 
		
- Output:
	- beta_hat: n_genes length vector of estimated regulated effect sizes, where beta_j != 0 indicates j th gene is not selected in HuberLasso regression.  

&emsp; &emsp;

### 11. HuberLasso_SP(X, y, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True)
***Estimate selection probability using HuberLasso function solving by APGD.***

- Input:
	- X: expression levels of n_genes target genes (TGs).
	- y: expression levels of a transcription factor (TF).
	- alphas: the grid sets of alpha (in [0,1]) used to calculate selection probabilities of genes.
	- n_lambda: the number of lambdas
	- B: the number of half-sample resampling used to calculate selection probabilities of genes.
		- default: 500
	- gamma: initial value of gamma in APGD.
		- default: 1000
	- niter: the maximum number of APGD to solve Huber Lasso regression.
		- default: 2000
	- crit_beta: converge criterion of change of beta.
		- default: 1e-4
	- crit_obj: converge criterion of change of objective function.
		- default: 1e-8
	- timer: decide if exist the output report.
		- default: True  

- Output:
	- sp_huberlasso: n_genes length vector of selection probability.

&emsp; &emsp;

### 12. HuberEnet_Beta(X, y, lambda0, alpha0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False)
***Estimate beta_hat using Huber Elastic Net function***

- Input:
	- X: expression levels of n_genes target genes (TGs).
	- y: expression levels of a transcription factor (TF).
	- lambda0: one of parameters in Huber Elastic Net regression, which controls the number of nonzero coefficients.
	- alpha0: one of parameters in Huber Elastic Net regression, which controls the numerical values of nonzero coefficients.
 	- method: the current methods must be 'APGD' or 'CVX'.
	- gamma: initial value of gamma in APGD.
		- default: 1000
	- niter: the maximum number of APGD to solve Huber Elastic Net regression.
		- default: 2000
 	- crit_beta: converge criterion of change of beta.
 	 	- default: 1e-4
 	- crit_obj: converge criterion of change of objective function.
		- default: 1e-8
	- quiet: decide if exist the output report.
		- default: FALSE  
		
- Output:
	- beta_hat: n_genes length vector of estimated regulated effect sizes, where beta_j != 0 indicates j th gene is not selected in Huber Elastic Net regression.

&emsp; &emsp;

### 13. HuberEnet_SP(X, y, alphas, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True)
***Estimate selection probability using HuberENET function solving by APGD.***

- Input:
	- X: expression levels of n_genes target genes (TGs).
	- y: expression levels of a transcription factor (TF).
	- alphas: the grid sets of alpha (in [0,1]) used to calculate selection probabilities of genes.
	- n_lambda: the number of lambdas
	- B: the number of half-sample resampling used to calculate selection probabilities of genes.
		- default: 500
	- gamma: initial value of gamma in APGD.
		- default: 1000
	- niter: the maximum number of APGD to solve Huber Elastic Net regression.
		- default: 2000
	- crit_beta: converge criterion of change of beta.
		- default: 1e-4
	- crit_obj: converge criterion of change of objective function.
		- default: 1e-8
	- timer: decide if exist the output report.
		- default: True  
		
- Output:
	- sp_huberenet: n_genes length vector of selection probability.  
 
&emsp; &emsp;

### 14. Enet_Beta(X, y, lambda0, alpha0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False)
***Estimate beta_hat using Elastic Net function.***

- Input:
	- X: expression levels of n_genes target genes (TGs).
	- y: expression levels of a transcription factor (TF).
	- lambda0: one of parameters in Elastic Net regression, which controls the number of nonzero coefficients.
	- alpha0: one of parameters in Elastic Net regression, which controls the numerical values of nonzero coefficients.
 	- method: the current methods must be 'APGD' or 'CVX'.
	- gamma: initial value of gamma in APGD.
		- default: 1000
	- niter: the maximum number of APGD to solve Elastic Net regression.
		- default: 2000
 	- crit_beta: converge criterion of change of beta.
 	 	- default: 1e-4
 	- crit_obj: converge criterion of change of objective function.
		- default: 1e-8
	- quiet: decide if exist the output report.
		- default: FALSE  
		
- Output:
	- beta_hat: n_genes length vector of estimated regulated effect sizes, where beta_j != 0 indicates j th gene is not selected in Elastic Net regression.

&emsp; &emsp;

### 15. Enet_SP(X, y, alphas, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True)
***Estimate selection probability using ENET function solving by APGD.***

- Input:
	- X: expression levels of n_genes target genes (TGs).
	- y: expression levels of a transcription factor (TF).
	- alphas: the grid sets of alpha (in [0,1]) used to calculate selection probabilities of genes.
	- n_lambda: the number of lambdas
	- B: the number of half-sample resampling used to calculate selection probabilities of genes.
		 - default: 500
	- gamma: initial value of gamma in APGD.
		- default: 1000
	- niter: the maximum number of APGD to solve Elastic Net regression.
		- default: 2000
	- crit_beta: converge criterion of change of beta.
		- default: 1e-4
	- crit_obj: converge criterion of change of objective function.
		- default: 1e-8
	- timer: decide if exist the output report.
		- default: True  
		
- Output:
	- sp_enet: n_genes length vector of selection probability.

&emsp; &emsp;

### 16. Lasso_Beta(X, y, lambda0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False)
***Estimate beta_hat using Lasso function.***

- Input:
	- X: expression levels of n_genes target genes (TGs).
	- y: expression levels of a transcription factor (TF).
	- lambda0: one of parameters in HuberNet regression, which controls the number of nonzero coefficients.
	- gamma: initial value of gamma in APGD.
		- default: 1000
	- niter: the maximum number of APGD to solve Lasso regression.
		- default: 2000
 	- crit_beta: converge criterion of change of beta.
 	 	- default: 1e-4
 	- crit_obj: converge criterion of change of objective function.
		- default: 1e-8
	- quiet: decide if exist the output report.
		- default: FALSE  
		
- Output:
	- beta_hat: n_genes length vector of estimated regulated effect sizes, where beta_j != 0 indicates j th gene is not selected in Lasso regression.

&emsp; &emsp;

### 17. Lasso_SP(X, y, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True)
***Estimate selection probability using Lasso function solving by APGD.***

- Input:
	- X: expression levels of n_genes target genes (TGs).
	- y: expression levels of a transcription factor (TF).
	- alphas: the grid sets of alpha (in [0,1]) used to calculate selection probabilities of genes.
	- n_lambda: the number of lambdas
	- B: the number of half-sample resampling used to calculate selection probabilities of genes.
		- default: 500
	- gamma: initial value of gamma in APGD.
		- default: 1000
	- niter: the maximum number of APGD to solve HuberNet regression.
		- default: 2000
	- crit_beta: converge criterion of change of beta.
		- default: 1e-4
	- crit_obj: converge criterion of change of objective function.
		- default: 1e-8
	- timer: decide if exist the output report.
		- default: True  

- Output:
	- sp_lasso: n_genes length vector of selection probability.

&emsp; &emsp;

### 18. Net_Beta(X, y, adj, lambda0, alpha0, method="APGD", gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, quiet=False)
***Estimate beta_hat using Net function.***

- Input:
	- X: expression levels of n_genes target genes (TGs).
	- y: expression levels of a transcription factor (TF).
	- adj: the adjacency matrix of network structure.
	- lambda0: one of parameters in Net regression, which controls the number of nonzero coefficients.
	- alpha0: one of parameters in Net regression, which controls the numerical values of nonzero coefficients.
 	- method: the current methods must be 'APGD' or 'CVX'.
	- gamma: initial value of gamma in APGD.
		- default: 1000
	- niter: the maximum number of APGD to solve HuberNet regression.
		- default: 2000
 	- crit_beta: converge criterion of change of beta.
 	 	- default: 1e-4
 	- crit_obj: converge criterion of change of objective function.
		- default: 1e-8
	- quiet: decide if exist the output report.
		- default: FALSE  
		
- Output:
	- beta_hat: n_genes length vector of estimated regulated effect sizes, where beta_j != 0 indicates j th gene is not selected in Net regression.

&emsp; &emsp;

### 19. Net_SP(X, y, adj, alphas, n_lambda, B=500, gamma=1000, niter=2000, crit_beta=1e-4, crit_obj=1e-8, timer=True)
***Estimate selection probability using Net function solving by APGD.***

- Input:
	- X: expression levels of n_genes target genes (TGs).
	- y: expression levels of a transcription factor (TF).
	- alphas: the grid sets of alpha (in [0,1]) used to calculate selection probabilities of genes.
	- n_lambda: the number of lambdas
	- B: the number of half-sample resampling used to calculate selection probabilities of genes.
		- default: 500
	- gamma: initial value of gamma in APGD.
		- default: 1000
	- niter: the maximum number of APGD to solve Net regression.
		- default: 2000
	- crit_beta: converge criterion of change of beta.
		- default: 1e-4
	- crit_obj: converge criterion of change of objective function.
		- default: 1e-8
	- timer: decide if exist the output report.
		- default: True 
		
- Output:
	- sp_net: n_genes length vector of selection probability.

&emsp; &emsp;

## Support functions:

### 1. Huber_gradient(delta, m)
***Calculate the gradient of Huber function for an input value delta.***

- Input:
  - delta: Input value delta.
  - m: Shape parameter, which is defaulted to be one-tenth of the interquartile range (IRQ).  
  
- Output:
  - value: The gradient of Huber function for an input value z. 

&emsp; &emsp;




