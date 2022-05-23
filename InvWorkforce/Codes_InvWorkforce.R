###############################
### Training Set Generation ###
############################### 

# Required packages
library(lhs)			# For Latin hypercube sampling
library(randomForest)	# For metamodel training

# In this example, we use a modifided version of pysd2r package instead of the original one.
source("path/to/ipysd.R")

# Connect Vensim to R
# You may need to run the following line twice.
py <- pysd_connect()
py <- read_vensim(py, "path/to/InventoryWorkforce.mdl")
py		# Check the status

# RMSE (Root Mean Square Error) function
rmse <- function(real, pred){
	error <- pred - real
	return(sqrt(mean(error^2)))
}

# Generate the training set
n <- 200
set.seed(111)
dat <- maximinLHS(n = n, k = 11, dup = 5)
tr <- matrix(NA, nrow = nrow(dat), ncol = (ncol(dat) + 1))
tr[, 1] <- qunif(dat[, 1], 50, 150)
tr[, 2] <- qunif(dat[, 2], 4, 12)		
tr[, 3] <- qunif(dat[, 3], 6, 18)
tr[, 4] <- qunif(dat[, 4], 9.5, 28.5)
tr[, 5] <- qunif(dat[, 5], 4, 12)
tr[, 6] <- qunif(dat[, 6], 1, 3)
tr[, 7] <- qunif(dat[, 7], 0.125, 0.375)
tr[, 8] <- qunif(dat[, 8], 1, 3)
tr[, 9] <- qunif(dat[, 9], 20, 60)
tr[, 10] <- qunif(dat[, 10], 2, 6)
tr[, 11] <- qunif(dat[, 11], 3, 9)
colnames(tr) = c("AvgDurEmp", "AvgTimeFillVac", "InvAdjTime", "LabAdjTime", "ManCycTime",
				"MinOrdProTime", "Productivity", "SafStoCov", "StdWorkWeek", "VacAdjTime", "WIPAdjTime",
				"penalty")
tr <- round(tr, 4) 

inv_tr <- NULL
desinv_tr <- NULL
for (i in 1:nrow(tr)){

	set_components(py, list("AvgDurEmp" = tr[i, 1], "AvgTimeFillVac" = tr[i, 2], "InvAdjTime" = tr[i, 3], "LabAdjTime" = tr[i, 4],
							"ManCycTime" = tr[i, 5], "MinOrdProTime" = tr[i, 6], "Productivity" = tr[i, 7], "SafStoCov" = tr[i, 8], "StdWorkWeek" = tr[i, 9],
							"VacAdjTime" = tr[i, 10], "WIPAdjTime" = tr[i, 11]))
	out <- run_model_v2(py, vals = list("Inv" = 10000 * (tr[i, 8] + tr[i, 6]),
										"WIPInv" = 10000 * tr[i, 5], 
										"Labor" = 10000 / (tr[i, 7] * tr[i, 9]),
										"Vac" = tr[i, 2] * (10000 / (tr[i, 7] * tr[i, 9])) / tr[i, 1]))
	inv_tr <- rbind(inv_tr, out$Inv)
	desinv_tr <- rbind(desinv_tr, out$DesInv)
}

for (i in 1:nrow(inv_tr)){
	tr[i, 12] <- rmse(inv_tr[i, ], desinv_tr[i, ])
} 

tr <- as.data.frame(tr)
tr$penalty <- log(tr$penalty)

###########################
### Test Set Generation ###
###########################

set.seed(222)
test <- matrix(NA, nrow = 2000, ncol = 12)
test[, 1] <- runif(1000, 50, 150)
test[, 2] <- runif(1000, 4, 12)		
test[, 3] <- runif(1000, 6, 18)
test[, 4] <- runif(1000, 9.5, 28.5)
test[, 5] <- runif(1000, 4, 12)
test[, 6] <- runif(1000, 1, 3)
test[, 7] <- runif(1000, 0.125, 0.375)
test[, 8] <- runif(1000, 1, 3)
test[, 9] <- runif(1000, 20, 60)
test[, 10] <- runif(1000, 2, 6)
test[, 11] <- runif(1000, 3, 9)
colnames(test) = c("AvgDurEmp", "AvgTimeFillVac", "InvAdjTime", "LabAdjTime", "ManCycTime",
				"MinOrdProTime", "Productivity", "SafStoCov", "StdWorkWeek", "VacAdjTime", "WIPAdjTime",
				"penalty")
test <- round(test, 4) 

inv_test <- NULL
desinv_test <- NULL
for (i in 1:nrow(test)){

	set_components(py, list("AvgDurEmp" = test[i, 1], "AvgTimeFillVac" = test[i, 2], "InvAdjTime" = test[i, 3], "LabAdjTime" = test[i, 4],
							"ManCycTime" = test[i, 5], "MinOrdProTime" = test[i, 6], "Productivity" = test[i, 7], "SafStoCov" = test[i, 8], "StdWorkWeek" = test[i, 9],
							"VacAdjTime" = test[i, 10], "WIPAdjTime" = test[i, 11]))
	out <- run_model_v2(py, vals = list("Inv" = 10000 * (test[i, 8] + test[i, 6]),
										"WIPInv" = 10000 * test[i, 5], 
										"Labor" = 10000 / (test[i, 7] * test[i, 9]),
										"Vac" = test[i, 2] * (10000 / (test[i, 7] * test[i, 9])) / test[i, 1]))
	inv_test <- rbind(inv_test, out$Inv)
	desinv_test <- rbind(desinv_test, out$DesInv)
}

for (i in 1:nrow(inv_test)){
	test[i, 12] <- rmse(inv_test[i, ], desinv_test[i, ])
} 

test <- as.data.frame(test)
test$penalty <- log(test$penalty)

##################################################
### Metamodeling for a Continuous Model Output ###
##################################################

# If you would like to switch to the metamodeling step directly, you can use the provided training and test sets.
tr <- read.csv("path/to/tr_InvWorkforce.csv")
test <- read.csv("path/to/test_InvWorkforce.csv")

# Train the RF metamodel 
rf <- randomForest(penalty ~., data = tr, ntree = 100)

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

pred <- predict(rf, newdata = tr[, 1:11], nodes = TRUE, predict.all = TRUE)
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
	colnames(binarym) <- which(g$status == -1)
	for(j in 1:dim(tr)[1]){
		binarym[j, as.character(nodeids[j, i])] <- 1
	}
	for(k in 1:nofrules){
		cost <- c(cost, rmse(tr[which(binarym[, k] == 1), "penalty"], individual[which(binarym[, k] == 1), i]))
	}
	binarym <- replace(binarym, which(binarym == -1), 1)
	colnames(g) <- sapply(colnames(g), collapse)
	condset <- append(condset, getRules(g))
	finalbinarym <- cbind(finalbinarym, binarym)
}

# Solve the set partitioning problem to extract rules

library(gurobi)

model <- list()
model$A          <- finalbinarym
model$obj		 <- 1 + cost / max(cost)
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
