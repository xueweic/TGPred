<!-- README.md is generated from README.Rmd. Please edit that file -->

# APGD v.0.1.0

<!-- badges: start -->

The Accelerated Proximal Gradient Descent (APGD) algorithm is a R package to solve the penalized regression models, including 

- **HuberNet**: Huber loss function along with Network-based penalty function;
- **HuberLasso**: Huber loss function along with Lasso penalty function;
- **HuberENET**: Huber loss function along with Elastic Net penalty function;
- **ENET**: Mean square error loss function along with Elastic Net penalty function;
- **Lasso**: Mean square error loss function along with Lasso penalty function;
- **Net**: Mean square error loss function along with Network-based penalty function.

## Functions

<!-- badges: start -->

### 1. read_file(file_path)
**Read file path to get data**

- input:   
	- file_path: the path of the file (.txt .csv) and separator by tab ('\t'). 
  
- output:   
 - df: data frame    

&emsp; &emsp;
   	
### 2. ConstructNetwork(n_genes, structure)
**Construct the network structure from either Hierarchical Network or Barabasi-Albert Network in simulation studies**   

- input:   
 - n_genes: the number of genes   
 - structure: "HN": Hierarchical Network or  "BAN": Barabasi-Albert Network   

- output:   
 - adj_all: n_genes * n_genes dimensional symmetric adjacency matrix of network structure.	  
	

	
	
	
