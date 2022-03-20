### Lab 1: Causal Inference

## Part 1: Data and packages
library(tidyverse)

install.packages("CondIndTests")
install.packages("dHSIC")
install.packages("ppcor")

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

BiocManager::install("graph")
BiocManager::install("RBGL")
BiocManager::install("Rgraphviz")

install.packages("pcalg")

## Part 2: Graphs and Adjacency Matrices

library(qgraph)

## Collider-structure with qgraph
varnames <- c("X","Y","Z")
Adj <- matrix(c(0,0,1,
                0,0,1,
                0,0,0), 3,3, byrow = TRUE,
              dimnames = list(varnames,varnames))

qgraph(Adj, 
       labels = c("X","Y","Z"), # not necessary if Adj has dimnames
       #you can provide a custom layout by giving the x-y co-ordinates of each node
       layout = rbind(c(-1,-1),
                      c(1,-1),
                      c(0,1)))

## Fork-structure with qgraph
Fork_Adj <- matrix(c(0,0,0,
                     0,0,0,
                     1,1,0), 3,3, byrow = TRUE,
                   dimnames = list(varnames, varnames))

qgraph(Fork_Adj,
       layout = rbind(c(-1,-1),
                      c(1,-1),
                      c(0,1)))

## Mediator-structure with qgraph
Mediator_Adj <- matrix(c(0,0,1,
                         0,0,0,
                         0,1,0), 3,3, byrow = TRUE,
                       dimnames = list(varnames, varnames))

qgraph(Mediator_Adj,
       layout = rbind(c(-1,-1),
                      c(1,-1),
                      c(0,1)))

library(ggdag)

## Collider-structure with ggdag
coldag <- dagify(
  Z ~ X + Y,
  exposure = "X", # the "cause" variable you are interested in
  outcome = "Y", # the "effect" variable you are interested in
  # optional: give co-ordinates of the variables in the plot
  coords = list(x = c(X = -1, Y = 1, Z = 0),
                y = c(X = 0, Y = 0, Z = 1)) 
) 
ggdag_status(coldag) + theme_dag()

## Fork-structure with ggdag
forkdag <- dagify(
  X ~ Z,
  Y ~ Z,
  exposure = "X",
  outcome = "Y",
  coords = list(x = c(X = -1, Y = 1, Z = 0),
                y = c(X = 0, Y = 0, Z = 1))
)

ggdag_status(forkdag) + theme_dag()

## Mediator-structure with ggdag
meddag <- dagify(
  Y ~ X,
  Z ~ Y,
  exposure = "X",
  outcome = "Z",
  coords = list(x = c(X = -1, Y = 1, Z = 3),
                y = c(X = 0, Y = 0, Z = 0))
)

ggdag_status(meddag) + theme_dag()

## D-separation and DAGs

eddag <- dagify(
  EA ~ CI,
  AI ~ CI + EA + U,
  Inc ~ EA + CI + AI + U ,
  exposure = "EA", # the "cause" variable you are interested in
  outcome = "Inc", # the "effect" variable you are interested in
  # optional: give co-ordinates of the variables in the plot
  coords = list(x = c(EA = -1,CI = 0,AI =0 ,Inc =1 , U = -1),
                y = c(EA = 0,CI =1 ,AI = -1,Inc = 0, U = -1)) 
) 
ggdag_status(eddag) + theme_dag()

## All paths between EA and IE
    ## EA --> Inc
    ## EA --> AI --> Inc
    ## EA --> CI --> Inc
    ## EA --> CI --> AI --> Inc
    ## EA --> AI --> CI --> Inc
    ## EA --> AI --> U --> Inc
    ## EA --> CI --> AI --> U --> Inc

## If you want to estimate the causal effect of EA on Inc, you should block CI
## because this is both a backdoor path and a confounder. 
## The paths that transmit causal associations are the direct ones and the ones through AI.
## The path from EA --> AI --> U --> Inc is already blocked, because
## AI is a collider between EA and U.

## The researcher should condition on CI, see reasons above.

library(dagitty)
adjustmentSets(eddag)

## Measure EA, Inc, CI for sure. And if wanted you can measure AI
## if you're interested in the moderating effect of it and the indirect effect.

## Part 4: Predictive Analysis

mdata <- readRDS("mdata.RDS")

## Best prediction model:

summary(lm(mal ~ Income + health + net, data = mdata))

## Interpretation: net has a significant effect on predicting mal even after
## controlling for both Income and health. (But smaller effect than income)

maldag <- dagify(
  mal ~ Income + net,
  net ~ Income,
  health ~ Income + mal + net,
  exposure = "net",
  outcome = "mal"
)

ggdag_status(maldag) + theme_dag()

## In order to estimate the causal effect of net on mal: only control on income
## as this gives a backdoor path into net. Furthermore, health = collider, so
## we should not condition on this.

adjustmentSets(maldag)

summary(lm(mal ~ net + Income, data = mdata))

## Comparison to the previous model: the effect of income has become smaller
## and more variance is explained in mal. Furthermore, effect of net has
## become twice as large.

## Part 5: Valid Adjustment Sets

library(CondIndTests)
library(dHSIC)

data <- readRDS("ex5_data.RDS")

Adj <- rbind(c(0,0,0,1,0,0,0,0,0), 
             c(0,0,1,1,0,0,0,0,0), 
             c(0,0,0,0,0,0,0,1,0), 
             c(0,0,0,0,1,1,0,0,0), 
             c(0,0,0,0,0,0,0,0,0), 
             c(0,0,0,0,0,0,1,1,0), 
             c(0,0,0,0,0,0,0,0,0), 
             c(0,0,0,0,0,0,0,0,1), 
             c(0,0,0,0,0,0,0,0,0))

names <- c("C", "A", "K", "X", "F", "D", "G", "Y", "H")
dimnames(Adj) = list(names,names)

laymat <- matrix(
  c(-1,   1,
    -.5,   1,
    .5,   1,
    -.75,   0,
    -.75,  -1,
    0,   0,
    0,  -1,
    1,   0,
    1,  -1),9,2,byrow = T)

vsize =15; esize = 10; asize = 10

qgraph(Adj, 
       layout = laymat, 
       vsize =vsize, esize = esize, asize = asize)

## Null-hypothesis: C and A are marginally independent

dhsic.test(data[,"C"], data[,"A"])

## Fail to reject null hypothesis, so we can safely assume they are marginally independent

## Null-hypothesis: X and G are independent given D

CondIndTest(data[,"X"], data[,"G"], data[,"D"])

## Fail to reject null hypothesis that they are independent!

## Two true conditional independence statements:
      ## A and F are independent given X
CondIndTest(data[, "A"], data[, "F"], data[, "X"]) ## true!
      ## G and Y are independent given D
CondIndTest(data[, "G"], data[, "Y"], data[, "D"]) ## true!

## Two false conditional independence statements:
      ## G and Y are independent
dhsic.test(data[,"G"], data[,"Y"]) ## not true! there is a spurious association when not controlling for D
      ## A and F are independent
dhsic.test(data[,"A"], data[,"F"]) ## not true! there is an association between A and F that is explained by X

## Finding two valid adjustment sets of X on Y:
    ## It is not necessary to condition on C (but you can do it), as this only has a relationship with X
    ## However, from A (or K) there is a backdoorpath into X. So we must condition either on A or K
    ## in any of our adjustment sets

summary(lm(Y ~ X + A, data = as.data.frame(data)))
summary(lm(Y ~ X + K, data = as.data.frame(data)))

## The true effect is both times around 1.94! So this is correct (true value = 2)

## Getting the dag through dagify

datadag <- dagify(
  H ~ Y,
  F ~ X,
  G ~ D,
  Y ~ D + K,
  K ~ A,
  X ~ C + A,
  D ~ X,
  exposure = "X",
  outcome = "Y",
  coords = list(x = c(X = 0, F = 0, D = 1, G = 1, Y = 3, H = 3, A = 0.5, C = -0.5, K = 2),
                y = c(X = 1, F = 0, D = 1, G = 0, Y = 1, H = 0, A = 2, C = 2, K = 2))
)

ggdag_status(datadag) + theme_dag()

adjustmentSets(datadag, type = "all") ## Mine were correct!

## Part 6: Generating data from an SCM

## Creating the DAG
scmdag <- dagify(
  X ~ Z,
  Y ~ X + Z,
  exposure = "X",
  outcome = "Y",
  coords = list(x = c(X = 0, Y = .5, Z = 1),
                y = c(X = 0, Z = 0, Y = .5))
)

ggdag_status(scmdag) + theme_dag()

## Generating data

n <- 1000

set.seed(3)

Z <- rnorm(n = n, mean = 0, sd = 1)
X <- 2*Z + rnorm(n = n, mean = 0, sd = 1)
Y <- X + 2*Z + rnorm(n = n, mean = 0, sd = 1)

## Estimating causal effect: Because Z has a fork structure, it is a confounder
## so to estimate the effect of X on Y, we should condition on Z

summary(lm(Y ~ X + Z))

## The effect is large and total explained variance is 95%

## Part 7: Causal Discovery using Conditional Independence Methods

data <- as.data.frame(readRDS("data_cd_ex1.RDS"))

## Discovering (in)dependencies

library(ppcor)

## Marginal densities
cor.test(data$X1, data$X2)
cor.test(data$X1, data$X3)
cor.test(data$X1, data$X4)
cor.test(data$X2, data$X3)
cor.test(data$X2, data$X4)
cor.test(data$X3, data$X4)

## All are dependent

## Conditional:

pcor.test(data$X1, data$X2, data$X3)
pcor.test(data$X1, data$X2, data$X4)
pcor.test(data$X1, data$X3, data$X2)
pcor.test(data$X1, data$X3, data$X4)
pcor.test(data$X1, data$X4, data$X2)
pcor.test(data$X1, data$X4, data$X3)
pcor.test(data$X2, data$X3, data$X1) ## INDEPENDENT
pcor.test(data$X2, data$X3, data$X4)
pcor.test(data$X2, data$X4, data$X1)
pcor.test(data$X2, data$X4, data$X3)
pcor.test(data$X3, data$X4, data$X1)
pcor.test(data$X3, data$X4, data$X2)

## all dependent except 1

pcor(data)

## Two independencies:
## X2 on X3, conditioned on X1
## X1 and X4, conditioned on X2 and X3

## If X1 has fork structure, gives different results than a mediator structure
## hence it makes sense that it implies the causal effects.

## Skeleton of the DAG:

adj_full <- matrix(1,4,4)
diag(adj_full) <- 0

# make the layout custom (optional)
layout = matrix(c(0,1,-1,0,1,0,0,-1),4,2,byrow = T)

# Make the ``full'' graph
qgraph(adj_full, labels = names, layout = layout, directed = FALSE, title = "Full Undirected Graph", title.cex = 1.25, vsize = 15)

# Remove the edges between X2 - X3 and X1- x4
adj_full <- matrix(1,4,4)
diag(adj_full) <- 0
adj <- adj_full
adj[2,3] <- adj[3,2] <- 0
adj[1,4] <- adj[4,1] <- 0

# make the layout custom (optional)
layout = matrix(c(0,1,-1,0,1,0,0,-1),4,2,byrow = T)

par(mfrow=c(1,2))
qgraph(adj_full, labels = names, layout = layout, directed = FALSE, title = "Full Undirected Graph", title.cex = 1.25, vsize = 15)
qgraph(adj, labels = names, layout = layout, directed = FALSE, title = "Estimated Skeleton", title.cex = 1.25, vsize = 15)

## CP-DAG

## There are four triplets you must consider:  
##  A) X2 - X1 - X3  
##  B) X1 - X3 - X4  
##  C) X1 - X2 - X4  
##  D) X2 - X4 - X3

## We rule out A) because we found above that  
##  X2 and X3 are independent given X1

## We rule out C) and B) because we found that   
##  X1 is independent of X4 given {X2 , X3}

## But X2 and X3 are always dependent (given either X4 or {X1, X4})  
##  That means X2 -> X4 <- X3

cpdag <- adj
cpdag[4,2] <- 0 # we know that this arrow goes X2 -> X4
cpdag[4,3] <- 0 # we know the direction is X3 -> X4

# extra touch - making a matrix that indicates we have a mix of directed (true) and undirected (false) edges
cptf <- matrix(FALSE, 4,4)
cptf[2,4] <- cptf[3,4] <- TRUE

par(mfrow = c(1,1))
qgraph(cpdag, labels = names, layout = layout, directed = cptf, title = "Estimated CPDAG", title.cex = 1.25, asize = 8, vsize = 15)

## Markov Equivalenc Class

dag1 <- matrix(c(
  0  ,  0  ,  1  ,  0,
  1  ,  0  ,  0  ,  1,
  0  ,  0  ,  0  ,  1,
  0  ,  0  ,  0  ,  0), 4, 4, byrow = T)

dag2 <- matrix(c(
  0  ,  1  ,  0  ,  0,
  0  ,  0  ,  0  ,  1,
  1  ,  0  ,  0  ,  1,
  0  ,  0  ,  0  ,  0), 4, 4, byrow = T)

dag3 <- matrix(c(
  0  ,  1  ,  1  ,  0,
  0  ,  0  ,  0  ,  1,
  0  ,  0  ,  0  ,  1,
  0  ,  0  ,  0  ,  0), 4, 4, byrow = T)

par(mfrow = c(1,3))
qgraph(dag1, labels = names, layout = layout, directed = TRUE, asize = 8, vsize = 15)
qgraph(dag2, labels = names, layout = layout, directed = TRUE, title = "Estimated Markov Equiv. Class", title.cex = 1.25, asize = 8, vsize = 15)
qgraph(dag3, labels = names, layout = layout, directed = TRUE, asize = 8, vsize = 15)

## These DAGS in the markov equivalence class also imply different causal effects.


## 7.2: PC-Algorithm

library(pcalg)

suffStat <- list(C = cor(data), n = nrow(data))

pc_fit1 <- pc(suffStat = suffStat, indepTest = gaussCItest,
              p = ncol(data), alpha = 0.01)
# This is the default plotting method for pcalg - uses Rgraphviz
plot(pc_fit1, main = "Inferred CPDAG using pcalg")

## Now we have the CPDAG
## To get the Markov Equivalence Class:

# Extract the adjacency matrix of the cpdag from pc_fit1
cpdag_mat <- as(pc_fit1,"matrix")

# Each row is a DAG adjacency matrix in vector form (by rows)
res1 <- pdag2allDags(cpdag_mat)

# We can get the adjacency matrix of an individual DAG using
res1_dags <- list()
for(i in 1:nrow(res1$dags)){
  res1_dags[[i]] <- t(matrix(res1$dags[i,],4,4,byrow = TRUE))
}
# Notice we have to transpose the adjacency matrix here for qgraph!

# We can plot each of these just as we did above
par(mfrow = c(1,3))
for(i in 1:3){
  qgraph(res1_dags[[i]], labels = names, layout = layout, directed = TRUE, asize = 8, vsize = 15)
}

## ida() function to estimate causal effect according to each DAG in Markov Equivalence set

ida(1,4,cov(data), pc_fit1@graph, verbose = TRUE)
