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


## Merged dataset
# Merge the year data for both caseids in allcites
allcites_with_years <- allcites %>%
  left_join(judicial %>% select(caseid, year), by = c("V1" = "caseid")) %>%
  rename(year1 = year) %>%
  left_join(judicial %>% select(caseid, year), by = c("V2" = "caseid")) %>%
  rename(year2 = year)

head(allcites_with_years)

## Figure 6

# Create year range and initialize storage lists
year_range <- 1799:2002

# Fig 6 ids
fig6_ids <- c("25347", "21109")

subsets <- list()
authority_scores_by_year <- list()

# Process each year in reverse
for (yr in rev(year_range)) {
  # Subset data for years within the range
  subsets[[as.character(yr)]] <- subset(allcites_with_years, year1 >= 1799 & year1 <= yr)
  complete_subsets <- subsets[[as.character(yr)]]

  # Create graph for the full network
  nodes_complete <- unique(c(complete_subsets$V1, complete_subsets$V2))
  vertices <- data.frame(name = nodes_complete)
  g_complete <- graph_from_data_frame(complete_subsets, directed = TRUE, vertices = vertices)

  # Run HITS algorithm
  hits_result_complete <- authority_score(g_complete)

  # Authority scores
  authority_scores_complete <- hits_result_complete$vector

  # Normalize authority scores
  normalized_authority_complete <- authority_scores_complete / sqrt(sum(authority_scores_complete^2, na.rm = TRUE))

  # Store scores for the year
  authority_scores_by_year[[as.character(yr)]] <- normalized_authority_complete

  # Display specific case IDs' authority scores for the year
  if (all(fig6_ids %in% names(normalized_authority_complete))) {
    cat(sprintf("\nYear: %d\n", yr))
    print(normalized_authority_complete[fig6_ids])
  }
}
