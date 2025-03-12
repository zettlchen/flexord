library("flexmix")
library("flexord")
library("flexclust")

# example data
data("iris", package = "datasets")
my_iris <- subset(iris, select=setdiff(colnames(iris), "Species")) |>
    as.matrix()

# cluster one model with a scale parameter similar to the default for 3 components
m1 <- stepFlexmix(my_iris~1,
                 model=FLXMCregnorm(zeta_p=c(0.23, 0.06, 1.04, 0.19)),
                 k=3)

summary(m1)

# rand index of clusters vs species
randIndex(clusters(m1), iris$Species)

# cluster one model with default scale parameter
m2 <- stepFlexmix(my_iris~1,
                 model=FLXMCregnorm(G=3),
                 k=3)

summary(m2)

# rand index of clusters vs species
randIndex(clusters(m2), iris$Species)


# rand index between both models (should be around 0.8)
randIndex(clusters(m1), clusters(m2))
