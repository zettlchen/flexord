library("flexmix")
library("flexord")
library("flexclust")


set.seed(0xdeaf)

# Sample data
k <- 4     # nr of clusters
nvar <- 10  # nr of variables
size <- sample(1:6, size=nvar, replace=TRUE)  # nr of trials 
N <- 100   # obs. per cluster


# random probabilities per component
probs <- lapply(seq_len(k), \(ki) runif(nvar, 0.01, 0.99))

# sample data
dat <- lapply(probs, \(p) {
    mapply(\(p_i, size_i) {
        rbinom(N, size_i, p_i)
    }, p, size, SIMPLIFY=FALSE) |> do.call(cbind, args=_)
}) |> do.call(rbind, args=_)

true_clusters <- rep(1:4, rep(N, k))

# Sample data is drawn from a binomial distribution but we fit
# a multinomial meaning the model is mis-specified.
# Note that for the multinomial distribution we expect values to lie inside
# 1:(size+1) hence we add +1.

# Cluster without regularization
m1 <- stepFlexmix((dat+1L)~1, model=FLXMCregmultinom(r=size+1L, alpha=0), k=k)

# Cluster with regularization
m2 <- stepFlexmix((dat+1L)~1, model=FLXMCregmultinom(r=size+1L, alpha=1), k=k)

# Both models are mostly able to reconstruct the true clusters (ARI ~ 0.95)
# (it's a very easy clustering problem)
# Small values for the regularization don't seem to affect the ARI (much)
randIndex(clusters(m1), true_clusters)
randIndex(clusters(m2), true_clusters)
