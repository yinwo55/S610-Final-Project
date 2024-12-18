library(tidyverse)
library(microbenchmark)
library(dplyr)
library(stringr)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(gridExtra)
library(xtable)
library(Matrix)
library(foreach)
library(doParallel)
library(tidyr)
options(digits = 3)


### 1. Table 1
# Loading Data
allcite_path <- "C:\\Users\\yonwo\\OneDrive\\바탕 화면\\PH.D\\course work\\3rd Semester\\Intro to Statistical Computing\\Final Project\\Judicial\\allcites.txt"
judicial_path <- "C:\\Users\\yonwo\\OneDrive\\바탕 화면\\PH.D\\course work\\3rd Semester\\Intro to Statistical Computing\\Final Project\\Judicial\\judicial.csv"
indeg_path <- "C:\\Users\\yonwo\\OneDrive\\바탕 화면\\PH.D\\course work\\3rd Semester\\Intro to Statistical Computing\\Final Project\\Judicial\\indegmat.txt"
outdeg_path <- "C:\\Users\\yonwo\\OneDrive\\바탕 화면\\PH.D\\course work\\3rd Semester\\Intro to Statistical Computing\\Final Project\\Judicial\\outdegmat.txt"

allcite <- read.csv(allcite_path, header = FALSE, sep = " ")
judicial <- read.csv(judicial_path, header = TRUE, sep = ",")
indeg <- read.csv(indeg_path, header = TRUE, sep = ",")
outdeg <- read.csv(outdeg_path, header = TRUE, sep = ",")


##########################################################
##### Case 5(Abortion Landmark) with Matrix Calculation
##########################################################


##### Matrix calculation #####

# Define the nodes in the five-case network
five_case_nodes <- c(25347, 27633, 28354, 29003, 29459)

# Filter edges where both `V1` and `V2` belong to the five-case network
Roe_v_wade_five <- allcite %>% 
  filter(V1 %in% five_case_nodes, V2 %in% five_case_nodes)

# Define edges from the filtered data
edges <- data.frame(from = Roe_v_wade_five$V1, to = Roe_v_wade_five$V2)

# Define nodes (extract unique nodes from edges)
nodes <- unique(c(edges$from, edges$to))

# Create the graph
g <- graph_from_data_frame(edges, directed = TRUE, vertices = data.frame(name = nodes))

# Generate Adjacency matrix (dense matrix)
adj_matrix <- as_adjacency_matrix(g, sparse = FALSE)

# Transposed adjacency matrix
transposed_matrix <- t(adj_matrix)

# Matrix-based Calculation with Sparse Matrix Optimization
# Sparse adjacency matrix
adj_matrix_sparse <- as(adj_matrix, "dgCMatrix")  # Convert to sparse matrix format

# Transposed sparse matrix
transposed_matrix_sparse <- t(adj_matrix_sparse)

# A^T * A using sparse matrix
authority_matrix_sparse <- transposed_matrix_sparse %*% adj_matrix_sparse

# Eigen decomposition with sparse matrix
eig_result_sparse <- eigen(as.matrix(authority_matrix_sparse))

# Extracting the eigenvector corresponding to the largest eigenvalue
raw_authority_scores_sparse <- eig_result_sparse$vectors[, 1]

# Adjusting the sign
if (raw_authority_scores_sparse[1] < 0) {
  raw_authority_scores_sparse <- -raw_authority_scores_sparse
}

# Euclidean normalization
normalized_authority_scores_matrix <- raw_authority_scores_sparse / sqrt(sum(raw_authority_scores_sparse^2))

##### HITS calculation #####

# HITS algorithm using hits_scores()
hits_result <- hits_scores(g)

# Extract authority scores
authority_scores <- hits_result$authority

# Normalize authority scores (Euclidean normalization)
normalized_authority_scores_hits <- authority_scores / sqrt(sum(authority_scores^2))

# Result Comparison
results <- rbind(
  Matrix = normalized_authority_scores_matrix,
  HITS = normalized_authority_scores_hits
)

print(results)


# Benchmark the two methods
benchmark_results <- microbenchmark(
  Matrix_Calculation = {
    # Matrix-based calculation (optimized with sparse matrices)
    authority_matrix_sparse <- transposed_matrix_sparse %*% adj_matrix_sparse
    eig_result_sparse <- eigen(as.matrix(authority_matrix_sparse))
    raw_authority_scores_sparse <- eig_result_sparse$vectors[, 1]
  },
  HITS_Calculation = {
    # HITS algorithm
    hits_result <- hits_scores(g)
    authority_scores <- hits_result$authority
  },
  times = 100
)

# Print benchmark results
print(benchmark_results)


##########################################################
##### Authority score for Fig 6 by Using HITS Packages
##########################################################


# Extract all years
years <- colnames(indeg)[-1]  # Exclude the first column (caseid)

# Initialize a file for results
output_file <- "C:\\Users\\yonwo\\OneDrive\\바탕 화면\\PH.D\\course work\\3rd Semester\\Intro to Statistical Computing\\Final Project\\authority_hub_scores_per_year.txt"
write.csv(data.frame(), output_file, row.names = FALSE)  # Create an empty file

# Maximum size for submatrices
max_matrix_size <- 200  # Adjust based on memory constraints

# Function to process a submatrix
process_submatrix <- function(indegree, outdegree) {
  adj_matrix <- as(as.matrix(indegree) %*% t(as.matrix(outdegree)), "dgCMatrix")
  
  # Skip if adjacency matrix is empty
  if (sum(adj_matrix) == 0) {
    return(NULL)
  }
  
  # Create the graph from sparse matrix
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "directed", weighted = NULL)
  
  # Apply HITS algorithm
  hits_result <- hits_scores(g)
  
  # Extract authority and hub scores
  authority_scores <- hits_result$authority
  hub_scores <- hits_result$hub
  
  # Normalize the scores
  normalized_authority_scores <- authority_scores / sqrt(sum(authority_scores^2, na.rm = TRUE))
  normalized_hub_scores <- hub_scores / sqrt(sum(hub_scores^2, na.rm = TRUE))
  
  return(data.frame(
    caseid = rownames(indegree),
    authority_score = normalized_authority_scores,
    hub_score = normalized_hub_scores
  ))
}

# Function to process a single year
process_year <- function(year) {
  cat("Processing year:", year, "\n")
  
  # Extract indegree and outdegree for the selected year
  indegree <- indeg %>% select(caseid, !!sym(year))
  outdegree <- outdeg %>% select(caseid, !!sym(year))
  
  # Ensure numeric data and replace NAs with zeros
  indegree <- indegree %>% mutate(across(-caseid, ~ifelse(is.na(.), 0, as.numeric(.))))
  outdegree <- outdegree %>% mutate(across(-caseid, ~ifelse(is.na(.), 0, as.numeric(.))))
  
  # Filter rows/columns with non-zero values
  non_zero_rows <- which(rowSums(as.matrix(indegree[,-1])) > 0)
  non_zero_cols <- which(rowSums(as.matrix(outdegree[,-1])) > 0)
  common_indices <- intersect(non_zero_rows, non_zero_cols)
  
  if (length(common_indices) == 0) {
    warning(paste("Skipping year", year, "- no valid data."))
    return(NULL)
  }
  
  indegree <- indegree[common_indices, , drop = FALSE]
  outdegree <- outdegree[common_indices, , drop = FALSE]
  
  # Break the matrix into smaller chunks
  chunk_indices <- split(1:nrow(indegree), ceiling(seq_along(1:nrow(indegree)) / max_matrix_size))
  year_result <- data.frame()
  
  for (chunk in chunk_indices) {
    sub_indegree <- indegree[chunk, -1, drop = FALSE]
    sub_outdegree <- outdegree[chunk, -1, drop = FALSE]
    sub_result <- process_submatrix(sub_indegree, sub_outdegree)
    if (!is.null(sub_result)) {
      year_result <- bind_rows(year_result, sub_result)
    }
  }
  
  # Add year information
  year_result$year <- as.numeric(gsub("X", "", year))
  return(year_result)
}

# Process each year and save results to file
for (year in years) {
  year_result <- process_year(year)
  if (!is.null(year_result)) {
    write.table(year_result, output_file, sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
    gc()  # Force garbage collection to free memory
  }
}

cat("Processing complete. Results saved to", output_file, "\n")














