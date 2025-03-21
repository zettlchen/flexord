library("flexmix")
library("flexord")
library("flexclust")


set.seed(0xdeaf)

# Sample data
k <- 4     # nr of clusters
nvar <- 10  # nr of variables
r <- sample(2:7, size=nvar, replace=TRUE)  # nr of categories
N <- 100   # obs. per cluster


# random probabilities per component
probs <- lapply(seq_len(k), \(ki) runif(nvar, 0.01, 0.99))

# sample data by drawing from a binomial distribution with size = r - 1
# values are expect values to lie inside 1:r hence we add +1.
dat <- lapply(probs, \(p) {
    mapply(\(p_i, r_i) {
        rbinom(N, r_i, p_i) + 1
    }, p, r-1, SIMPLIFY=FALSE) |> do.call(cbind, args=_)
}) |> do.call(rbind, args=_)

true_clusters <- rep(1:4, rep(N, k))

# Cluster without regularization
m1 <- stepFlexmix(dat~1, model=FLXMCregmultinom(r=r, alpha=0), k=k)

# Cluster with regularization
m2 <- stepFlexmix(dat~1, model=FLXMCregmultinom(r=r, alpha=1), k=k)

# Both models are mostly able to reconstruct the true clusters (ARI ~ 0.95)
# (it's a very easy clustering problem)
# Small values for the regularization don't seem to affect the ARI (much)
randIndex(clusters(m1), true_clusters)
randIndex(clusters(m2), true_clusters)
