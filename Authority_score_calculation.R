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

allcites <- read.csv(allcite_path, header = FALSE, sep = " ")
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
Roe_v_wade_five <- allcites %>% 
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


# merging year with allcites
allcites_year <- merge(allcites, judicial[, c("caseid", "year")], 
                 by.x = "V1", by.y = "caseid", all.x = TRUE)

allcites_year <- allcites_year[, c("year", "V1", "V2")]

head(allcites_year)


# allcite drop_year using for loop
unique_years <- sort(unique(allcites_year$year), decreasing = TRUE)

year_range <- 1799:2002

subsets <- list()
results <- list()

for (yr in rev(year_range)) {
  subsets[[as.character(yr)]] <- subset(allcites_year, year >= 1799 & year <= yr)
  current_subset <- subsets[[as.character(yr)]]

  nodes_complete <- unique(c(current_subset$V1, current_subset$V2))
  vertices <- data.frame(name = nodes_complete)

  g_complete <- graph_from_data_frame(current_subset, directed = TRUE, vertices = vertices)
}

  edge_nodes <- unique(c(current_subset$V1, current_subset$V2))
  vertex_nodes <- vertices$name
  
  missing_nodes <- setdiff(edge_nodes, vertex_nodes)
  if (length(missing_nodes) > 0) {
    print(missing_nodes)
  }


# Error messages occurred:
# graph_from_data_frame(current_subset, directed = TRUE, vertices = vertices)에서 다음과 같은 에러가 발생했습니다: 
#  Some vertex names in edge list are not listed in vertex data frame


# setdiff(unique(c(current_subset$V1, current_subset$V2)), nodes_complete)
# print(vertices)
# print(nodes_complete)
# str(vertices$name)









