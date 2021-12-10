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

- input:   
	- file_path: the path of the file (.txt .csv) and separator by tab ('\t'). 
  
- output:   	
	- df: data frame    

&emsp; &emsp;
   	
### 2. ConstructNetwork(n_genes, structure)
**Construct the network structure from either Hierarchical Network or Barabasi-Albert Network in simulation studies.**   

- input:
	- n_genes: the number of genes   
	- structure: "HN": Hierarchical Network or  "BAN": Barabasi-Albert Network   

- output:   
	- adj_all: n_genes * n_genes dimensional symmetric adjacency matrix of network structure.	  
	
&emsp; &emsp;

### 3.GraphicalModel(adj, a1=-0.7, a2=-0.1, b1=0.1, b2=0.7)
**Simulate a covariance matrix from the specific graph (Adj) based on a Gaussian graphical model.**  

- Input
	- adj: the adjacency matrix of network structure.
	- a1, a2, b1, b2: parameters for constructing domain [a1, a2] union [b1, b2].  
	  default: a1 = -0.7, a2 = -0.1, b1 = 0.1, b2 = 0.7
- Output
	- sigma: covariance matrix of target genes according to network structure.






	
	
	
