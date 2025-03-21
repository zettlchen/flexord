library("flexmix")
library("flexord")
library("flexclust")

# example data
data("iris", package = "datasets")
my_iris <- subset(iris, select=setdiff(colnames(iris), "Species")) |>
    as.matrix()

# cluster one model with a scale parameter similar to the default for 3 components
params <- FLXMCregnorm_defaults(my_iris, zeta_p = c(0.23, 0.06, 1.04, 0.19))
m1 <- stepFlexmix(my_iris ~ 1, k = 3, 
    model=FLXMCregnorm(params = params))
summary(m1)

# rand index of clusters vs species
randIndex(clusters(m1), iris$Species)

# cluster one model with default scale parameter
params <- FLXMCregnorm_defaults(my_iris, k = 3)
m2 <- stepFlexmix(my_iris ~ 1, k = 3,
    model = FLXMCregnorm(params = params))
summary(m2)

# rand index of clusters vs species
randIndex(clusters(m2), iris$Species)

# rand index between both models (should be >= 0.8)
randIndex(clusters(m1), clusters(m2))
