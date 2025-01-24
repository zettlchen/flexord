\donttest{
library("flexmix")
library("flexord")
library("flexclust")


# Sample data
k = 4     # nr of clusters
size = 4  # nr of trials
N = 100   # obs. per cluster

set.seed(0xdeaf)

# random probabilities per component
probs = lapply(seq_len(k), \(ki) runif(10, 0.01, 0.99))

# sample data
dat = lapply(probs, \(p) {
    lapply(p, \(p_i) {
        rbinom(N, size, p_i)
    }) |> do.call(cbind, args=_)
}) |> do.call(rbind, args=_)

true_clusters = rep(1:4, rep(N, k))

# Sample data is drawn from a binomial distribution but we fit
# beta-binomial which is a slight mis-specification but the
# beta-binomial can be seen as a generalized binomial.

# Cluster without regularization
m1 = stepFlexmix(dat~1, model=FLXMCbetabinom(size=size, alpha2=0), k=k)

# Cluster with regularization
m2 = stepFlexmix(dat~1, model=FLXMCbetabinom(size=size, alpha2=1), k=k)

# Both models are mostly able to reconstruct the true clusters (ARI ~ 0.95)
# (it's a very easy clustering problem)
# Small values for the regularization don't seem to affect the ARI (much)
randIndex(clusters(m1), true_clusters)
randIndex(clusters(m2), true_clusters)
}
