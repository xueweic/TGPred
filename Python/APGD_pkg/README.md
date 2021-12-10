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
		- default: a1 = -0.7, a2 = -0.1, b1 = 0.1, b2 = 0.7.  
		
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
	
	
