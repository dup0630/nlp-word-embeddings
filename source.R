
load("svd_saved.RData")
library(tidyverse)

# Load and normalize
cooccur <- read.csv("data/co_occur.csv", header = FALSE)
dictionary <- read.csv("data/dictionary.txt", header = FALSE)$V1
M_scaled <- apply(cooccur, 1:2, \(x) log(1 + x, base = exp(1)))
M_scaled[1:8, 1:8]

# SVD:
svd_M <- svd(M_scaled, nu = 100, nv = 100)

s_values <- svd_M$d
right_s_vectors <- svd_M$u
rownames(right_s_vectors) <- dictionary
left_s_vectors <- svd_M$v
rownames(left_s_vectors) <- dictionary

# Compute the approximation:
M_approx <- right_s_vectors %*% diag(s_values) %*% t(left_s_vectors)
M_approx[1:8, 1:8]

plot(s_values)

s_values[1:20]

# Create empty matrices to store the top words of each vector.
pos_top_words <- matrix(rep(0, 10*100), nrow = 10)
neg_top_words <- matrix(rep(0, 10*100), nrow = 10)
colnames(pos_top_words) <- paste("v", 1:100, sep = "")
colnames(neg_top_words) <- paste("v", 1:100, sep = "")

# Fill the matrix with a loop.
for(i in 1:100){
  # We sort each vector and extract the first 10 row names.
  pos_top_words[,i] <- names(sort(right_s_vectors[, i], decreasing = TRUE))[1:10]
  neg_top_words[,i] <- names(sort(right_s_vectors[, i], decreasing = FALSE))[1:10]
}

pos_top_words[, 1:6]
neg_top_words[, 1:6]

example_vectors <- c(10, 35, 48, 60, 97)
pos_top_words[, example_vectors]
neg_top_words[, example_vectors]


# Normalize the rows of the singular vector matrix
vectors <- t(apply(right_s_vectors, 1, function(x){x / norm(x, "2")}))
rownames(vectors) <- dictionary

V_1 <- vectors["woman",]
V_2 <- vectors["man",]
V <- V_1 - V_2

words <- c("boy", "girl", "brother", "sister", "king", "queen", "he", "she",
           "john", "mary", "wall", "tree")
indices <- which(dictionary %in% words)

# This step is necessary to ensure that both indices and words are in the same
# corresponding order
words <- rownames(vectors[indices, ])

projections_V <- vectors[indices, ] %*% V

proj_tibble <- tibble(words, projections_V)

proj_tibble %>% 
  ggplot(aes(x = words, y = projections_V)) + 
  geom_col(fill = "#009E73") + 
  labs(y = "Projections", x = "Words")

words <- c("math", "matrix", "history", "nurse", "doctor", "pilot", "teacher", 
           "engineer", "science", "arts", "literature", "bob", "alice")
indices <- which(dictionary %in% words)

# This step is necessary to ensure that both indices and words are in the same 
# corresponding order.
words <- rownames(vectors[indices, ])

projections_V <- vectors[indices, ] %*% V

proj_tibble <- tibble(words, projections_V)

proj_tibble %>% 
  ggplot(aes(x = words, y = projections_V)) + 
  geom_col(fill = "#009E73") + 
  labs(y = "Projections", x = "Words") + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))


d_montreal <- function(word){
  vectors[word, ] %*% vectors["montreal", ]
}

similar_montreal <- sapply(dictionary, d_montreal) %>% 
  sort(decreasing = TRUE) %>% 
  names()

similar_montreal[1:30]

solve_analogy <- function(row_vector, n_of_guesses = 5){
  v1 <- vectors[row_vector[1], ]
  v2 <- vectors[row_vector[2], ]
  v3 <- vectors[row_vector[3], ]
  target <- v2 - v1 + v3
  
  # Obtain the indices for the chosen words and remove them from possible guesses
  indices <- which(dictionary %in% c(row_vector[1], 
                                     row_vector[2], 
                                     row_vector[3])
                   )
  possible_vectors <- vectors[-indices, ]
  
  # Obtain the nearest word embeddings to the target:
  distances <- apply(possible_vectors, 1, \(x) x %*% target)
  result <- sort(distances, decreasing = TRUE) %>% names()
  
  result[1:n_of_guesses]
}

# Function for testing solve_analogy on rows of a dataset:
test <- function(data, n = 1){
  guesses <- apply(data, 1, solve_analogy, n_of_guesses = n)
  
  # For n>1, apply will return a matrix with a column output for each row input.
  # We want to turn these column outputs back to rows, and we use t().
  if(n>1){
    guesses <- t(guesses)
    colnames(guesses) <- paste("guess", 1:n, sep = "_")
  }
  
  cbind(data, guesses)
}


analogy_task <- read.csv("data/analogy_task.txt", header = FALSE, sep = " ")
test(analogy_task[1:8, ], n=3)

test_results <- test(analogy_task, n = 5)


# Filter the results to obtain the rows where all five guesses were wrong:
failed_tasks <- test_results %>% 
  filter(!Reduce(`|`, lapply(select(., -V4), \(x) x == V4)))

head(failed_tasks)
dim(failed_tasks)
dim(analogy_task)


failed_tasks[c(902, 945, 1285),]

failed_tasks[c(1021, 1024, 1026, 1029, 1033, 1040, 1044, 1048, 1053), ]
