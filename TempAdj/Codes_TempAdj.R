###############################
### Training Set Generation ###
###############################

# Required packages
library(lhs)			# For Latin hypercube sampling
library(pysd2r)			# For connecting Vensim to R
library(randomForest)	# For metamodel training

# Load SMTS to classify SD model outputs
setwd("path/to/SMTSClassifier")
source("predict.SMTS.r")

# Connect Vensim to R
# You may need to run the following line twice.
py <- pysd_connect()
py <- read_vensim(py, "path/to/TempAdj.mdl")
py		# Check the status

# Generate the training set
n <- 200
set.seed(1)
design <- maximinLHS(n = n, k = 2)
tr <- matrix(NA, ncol = 2, nrow = n)
tr[, 1] <- qunif(design[, 1], 1, 21)
tr[, 2] <- qunif(design[, 2], 1, 21)

art <- NULL
for (i in 1:nrow(tr)){
	set_components(py, list("tad" = tr[i, 1], "pd" = tr[i, 2]))
	out <- run_model(py)
	art <- rbind(art, c(tr[i, 1], tr[i, 2], out$art))
}

pred_art <- predict.SMTS(newdata = as.data.frame(art[, 3:ncol(art)]), modelPath = "path/to/SMTSClassifier/model.rds")
tr <- cbind.data.frame(tr, pred_art$classPred)
colnames(tr) <- c("tad", "pd", "class")

###########################
### Test Set Generation ###
###########################

# Generate the test set
set.seed(1)
tad <- runif(1000, min = 1, max = 21)
pd <- runif(1000, min = 1, max = 21)
test <- cbind(tad, pd)

art <- NULL
for (i in 1:nrow(test)){
	set_components(py, list("tad" = test[i, 1], "pd" = test[i, 2]))
	out <- run_model(py)
	art <- rbind(art, c(test[i, 1], test[i, 2], out$art))
}

pred_art <- predict.SMTS(newdata = as.data.frame(art[, 3:ncol(art)]), modelPath = "path/to/SMTSClassifier/model.rds")
test <- cbind.data.frame(test, pred_art$classPred)
colnames(test) <- c("tad", "pd", "class")

#####################################
### Metamodeling for Class Labels ###
##################################### 

# If you would like to switch to the metamodeling step directly, you can use the provided training and test sets.
tr <- read.csv("path/to/tr_TempAdj.csv")
test <- read.csv("path/to/test_TempAdj.csv")

# Train the RF metamodel
rf <- randomForest(class ~ ., data = tr, ntree = 100) 

#######################
### Rule Extraction ###
#######################

# These codes are taken from the link given below:
# https://stats.stackexchange.com/questions/41443/how-to-actually-plot-a-sample-tree-from-randomforestgettree

##################################
### Return the rules of a tree ###
##################################

getRules <- function(tree){
	# store all rules into a list
	rules <- list()
	# start by the terminal nodes and find previous conditions
	id.leafs <- which(tree$status == -1)
	j <- 0
	for(i in id.leafs){
		j <- j + 1
		prevConds <- prevCond(tree,i)
		rules[[j]] <- prevConds$cond
		while(prevConds$id > 1){
			prevConds <- prevCond(tree, prevConds$id)
			rules[[j]] <- paste(rules[[j]], "&", prevConds$cond)
        }
		if(prevConds$id == 1){
			rules[[j]] <- paste(rules[[j]], "=>", tree$prediction[i])
		}
    }
	return(rules)
}

################################################
### Find the previous conditions in the tree ###
################################################

prevCond <- function(tree, i){
	if(i %in% tree$right_daughter){
		id <- which(tree$right_daughter == i)
		cond <- paste(tree$split_var[id], ">", tree$split_point[id])
	}
	if(i %in% tree$left_daughter){
		id <- which(tree$left_daughter == i)
		cond <- paste(tree$split_var[id], "<=", tree$split_point[id])
	}
	return(list(cond = cond, id = id))
}

###############################
### Remove spaces in a word ###
###############################

collapse <- function(x){
	x <- sub(" ", "_", x)
	return(x)
}

pred <- predict(rf, newdata = tr[, 1:2], nodes = TRUE, predict.all = TRUE)
nodeids <- matrix(NA, ncol = rf$ntree, nrow = dim(tr)[1])
individual <- matrix(NA, ncol = rf$ntree, nrow = dim(tr)[1])
nodeids <- attr(pred, "nodes")
individual <- pred$individual

finalbinarym <- NULL
condset <- NULL
cost <- NULL
for(i in 1:rf$ntree){
	g <- getTree(rf, k = i, labelVar = TRUE)
	nofrules <- length(which(g$status == -1))
	binarym <- matrix(0, nrow = dim(tr)[1], ncol = nofrules)
	colnames(binarym) = which(g$status == -1)
	for(j in 1:dim(tr)[1]){
		if (individual[j, i] == as.character(tr$class[j])){
			binarym[j, as.character(nodeids[j, i])] <- 1
		}
		else
			binarym[j, as.character(nodeids[j, i])] <- -1
	}
	for(k in 1:nofrules)
		cost <- cbind(cost, 1 - length(which(binarym[, k] == 1)) / length(which(binarym[, k] != 0)))
	binarym <- replace(binarym, which(binarym == -1), 1)
	colnames(g) <- sapply(colnames(g), collapse)
	condset <- append(condset, getRules(g))
	finalbinarym <- cbind(finalbinarym, binarym)
}

# Solve the set partitioning problem to extract rules

library(gurobi)
 
model <- list()
model$A          <- finalbinarym
model$obj		 <- 1 + cost
model$modelsense <- "min"
model$rhs        <- rep(1, dim(tr)[1])		
model$sense      <- rep("=", dim(tr)[1])	
model$vtype      <- rep("B", ncol(finalbinarym))		
params <- list(TimeLimit = 3600)

t1 <- proc.time()
result <- gurobi(model, params)
t2 <- proc.time()

print(result$objval)

# Extract rules and display them
extruleset <- condset[which(result$x >= 1 - 1e-05)]
extruleset
