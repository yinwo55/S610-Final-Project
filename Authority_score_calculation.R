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

# Combine authority scores into a data frame for plotting
authority_scores_df <- do.call(rbind, lapply(names(authority_scores_by_year), function(year) {
  scores <- authority_scores_by_year[[year]]
  data.frame(
    year = as.numeric(year),
    case_id = names(scores),
    authority_score = scores
  )
}))

# Filter for the specific case IDs of interest
filtered_df <- authority_scores_df[authority_scores_df$case_id %in% fig6_ids, ]

# Adjust the years after decision
years_since_decision <- list("25347" = 1973, "21109" = 1954)
filtered_df$years_after <- filtered_df$year - sapply(filtered_df$case_id, function(id) years_since_decision[[id]])

# Plot authority scores over time
ggplot(filtered_df, aes(x = years_after, y = authority_score, linetype = case_id)) +
  geom_line() +
  labs(
    title = "",
    x = "Years After Decision",
    y = "Authority Score"
  ) +
  annotate("text", x = 2, y = 0.06, label = "Roe v. Wade,\n310 U.S. 113 (1973)", fontface = "italic", hjust = 0) +
  annotate("text", x = 11, y = 0.02, label = "Brown v. Board of Education,\n347 U.S. 483 (1954)", fontface = "italic", hjust = 0) +
  scale_linetype_manual(values = c("25347" = "solid", "21109" = "dashed")) +
  scale_x_continuous(breaks = seq(0, 30, by = 5),
                     limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 0.065)) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(
    color = "black", fill = NA)
    )

## Figure 10

# Create year range and initialize storage lists
year_range <- 1799:2002

# Fig 10 ids
fig10_ids <- c("18501", "23115", "23601", "26918")

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
  if (all(fig10_ids %in% names(normalized_authority_complete))) {
    cat(sprintf("\nYear: %d\n", yr))
    print(normalized_authority_complete[fig10_ids])
  }
}

# Combine authority scores into a data frame for plotting
authority_scores_df <- do.call(rbind, lapply(names(authority_scores_by_year), function(year) {
  scores <- authority_scores_by_year[[year]]
  data.frame(
    year = as.numeric(year),
    case_id = names(scores),
    authority_score = scores
  )
}))

# Filter for the specific case IDs of interest
filtered_df_10 <- authority_scores_df[authority_scores_df$case_id %in% fig10_ids, ]

# Plot authority scores over time
ggplot(filtered_df_10, aes(x = year, y = authority_score, linetype = case_id)) +
  geom_line() +
  labs(
    title = "",
    x = "Year",
    y = "Authority Score"
  ) +
  annotate("text", x = 1945, y = 0.04, label = "Brown v. Mississippi,\n297 U.S. 278 (1936)",
           fontface = "italic", hjust = 0) +
  annotate("text", x = 1985, y = 0.05, label = "Miranda v. Arizona,\n384 U.S. 436 (1966)",
           fontface = "italic", hjust = 0) +
  annotate("text", x = 1973, y = 0.02, label = "Escobedo v. Illinois,\n378 U.S. 478 (1964)",
           fontface = "italic", hjust = 0) +
  annotate("text", x = 1984, y = 0.01, label = "Rhode Island v. Innis,\n446 U.S. 291 (1980)",
           fontface = "italic", hjust = 0) +
  labs(
    x = "Year",
    y = "Authority Score",
    title = "5th Amendment"
  ) +
  scale_linetype_manual(values = c("18501" = "dotted", "23115" = "dashed",
                                   "23601" = "solid", "26918" = "dotdash")) +
  scale_x_continuous(breaks = seq(1940, 2000, by = 10), limits = c(1940, 2000)) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(
    color = "black", fill = NA)
    )



##########################################################
######## Authority score testing by Using testthat
##########################################################

                                                     
## An Alternative Way to Calculate Hub and Authority Scores      
                                                     
# create case_year_map 
case_year_map <- setNames(judicial$year, judicial$caseid)

# define function: from case_year_map, filter only the data up to a specific year
subset_allcites_year <- function(allcites, yr) {
  V1_year <- case_year_map[as.character(allcites$V1)]
  valid_rows <- !is.na(V1_year) & V1_year <= yr
  return(allcites[valid_rows, ])
}

# Initialize year variables
start_year <- 1799
end_year <- 2002

# Initialize result storage lists
authority_scores_by_year <- list()
hub_scores_by_year <- list()

# Process data by year
for (yr in seq(end_year, start_year, by = -1)) {
  # Filter data by year
  allcites_subset <- subset_allcites_year(allcites, yr)
  
  # Define nodes
  nodes_complete <- unique(c(allcites_subset$V1, allcites_subset$V2))
  
  # Create the graph
  g_complete <- graph_from_data_frame(allcites_subset, directed = TRUE, vertices = data.frame(name = nodes_complete))
  
  # Run the HITS algorithm
  hits_result_complete <- hits_scores(g_complete)
  
  # Calculate Authority and Hub scores
  authority_scores_complete <- hits_result_complete$authority
  hub_scores_complete <- hits_result_complete$hub
  
  # Normalize scores
  normalized_authority_complete <- authority_scores_complete / sqrt(sum(authority_scores_complete^2))
  normalized_hub_complete <- hub_scores_complete / sqrt(sum(hub_scores_complete^2))
  
  # Store results by year
  authority_scores_by_year[[as.character(yr)]] <- normalized_authority_complete
  hub_scores_by_year[[as.character(yr)]] <- normalized_hub_complete
}

# Convert Authority scores to a data frame
auth_df <- data.frame()
hub_df <- data.frame()

for (yr in names(authority_scores_by_year)) {
  year_auth_scores <- data.frame(
    year = as.numeric(yr),
    caseid = names(authority_scores_by_year[[yr]]),
    authority_score = authority_scores_by_year[[yr]]
  )
  auth_df <- rbind(auth_df, year_auth_scores)

  year_hub_scores <- data.frame(
    year = as.numeric(yr),
    caseid = names(hub_scores_by_year[[yr]]),
    hub_score = hub_scores_by_year[[yr]]
  )
  hub_df <- rbind(hub_df, year_hub_scores)
}

# Output Authority scores for cases of interest                                                     
cases_of_interest_testing <- auth_df %>%
  filter(caseid %in% c(21109, 25347)) %>%
  mutate(
    years_after_decision = ifelse(caseid == 21109, year - 1954, year - 1973),
    case_name = case_when(
      caseid == 21109 ~ "Brown v. Board of Education, 347 U.S. 483 (1954)",
      caseid == 25347 ~ "Roe v. Wade, 310 U.S. 113 (1973)"
    )
  ) %>%
  filter(!is.na(authority_score))

cases_of_interest_testing <- cases_of_interest_testing %>%
  group_by(caseid) %>%
  mutate(years_after_decision = year - min(year)) %>%
  ungroup()

head(cases_of_interest_testing)



                                                     
## Testing 1
## Our computation with original data
library(testthat)
                                                     
# Filter df for the Figure 6 cases
filtered_df <- authority_scores_df[authority_scores_df$case_id %in% fig6_ids, ]

# Authmat subset to compare with Figure 6
caseids_test <- c(21109, 25347)
years_test <- which(colnames(authmat) == "X1954"):which(colnames(authmat) == "X2002")

authmat_test <- authmat[authmat[, 1] %in% caseids_test, c(1, years_test)]

print(authmat_test)

# Ensure years in authmat_test columns are properly formatted
colnames(authmat_test) <- as.numeric(gsub("^X", "", colnames(authmat_test)))

# Define a function to test matching authority scores
test_that("Authority scores match up to 3 decimal places, skipping NA values", {
  
  # Loop through each row of filtered_df to check authority scores
  for (i in 1:nrow(filtered_df)) {
    year <- filtered_df$year[i]
    case_id <- filtered_df$case_id[i]
    auth_score <- filtered_df$authority_score[i]
    
    # Skip rows with NA authority scores
    if (is.na(auth_score)) {
      next
    }
    
    # Find the matching column in authmat_test for the year
    year_col <- which(colnames(authmat_test) == as.character(year))
    
    # Check if year or case_id is invalid and skip
    if (length(year_col) == 0 || is.na(year_col)) {
      next
    }
    case_row <- which(rownames(authmat_test) == as.character(case_id))
    if (length(case_row) == 0 || is.na(case_row)) {
      next
    }
    
    # Check if the corresponding cell in authmat_test is NA and skip
    if (is.na(authmat_test[case_row, year_col])) {
      next
    }
    
    # Perform the test with tolerance for 3 decimal places
    expect_equal(
      authmat_test[case_row, year_col],
      auth_score,
      tolerance = 1e-3,
      info = paste("Mismatch at case_id:", case_id, "year:", year)
    )
  }
})


## Testing 2
## Comparison of Authority Scores Computed Using Different Approaches
                                                                                                   
# Define a function to test matching authority scores between two data frames
test_that("Authority scores match between filtered_df and cases_of_interest_testing", {
  
  # Loop through each row of filtered_df
  for (i in 1:nrow(filtered_df)) {
    year <- filtered_df$year[i]
    caseid <- filtered_df$case_id[i]
    auth_score <- filtered_df$authority_score[i]
    
    # Find the corresponding row in cases_of_interest_testing
    matching_row <- cases_of_interest_testing[
      cases_of_interest_testing$year == year & 
      cases_of_interest_testing$caseid == caseid, ]
    
    # Check if matching row exists and skip if not found
    if (nrow(matching_row) == 0) {
      next
    }
    
    # Perform the test with tolerance for 3 decimal places
    expect_equal(
      auth_score,
      matching_row$authority_score,
      tolerance = 1e-3,
      info = paste("Mismatch at caseid:", caseid, "year:", year)
    )
  }
})


                                                     
