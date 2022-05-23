tune.SMTS <- function(trainingdata, classes, tuningParamLevels = list(noftreelevels = c(10, 25, 50), nofnodelevels =c (10, 25, 50))){
	
	require(randomForest)
	require(plyr)
	
	# setwd("C:/Users/Mert/Desktop/R Package/R Package v2")
	
	if(is.data.frame(classes) == FALSE)
		stop("classes must be a data.frame!")
	if(is.data.frame(trainingdata) == FALSE)
		stop("trainingdata must be a data.frame!")
		
	trclass <- classes
	colnames(trclass) <- "x"
	uniqueclass <- unique(trclass$x)
	uniqueclass <- data.frame(x = uniqueclass, ID = c(1:length(uniqueclass)))
	trclass <- join(trclass, uniqueclass, type = "left")

	trainingdata <- as.matrix(cbind(trclass$ID, trainingdata))

	# source("trainingdata_preparation.r")
	
	#following codes obtains the training and test matrix in the paper
	#these are the ones that goes into random forest model, RFobs

	datatraintimestart <- proc.time()
	nofclass <- length(unique(trainingdata[, 1]))	# number of classes
	classtrain <- trainingdata[, 1] 				# classes of the training time series
	noftrain <- nrow(trainingdata) 					# number of training series
	seriesLen <- apply(trainingdata[, 2:ncol(trainingdata)], 1, function(x) sum(!is.na(x)))	# length of each series
	observations <- array(0, sum(seriesLen)-noftrain)	# observation array (storing all observations as a column)
	difference <- array(0, sum(seriesLen)-noftrain)		# difference array (storing difference between consecutive observations as a column)

	#for each time series observations should be standardized and
	#we need to concatenate the observation of each time series as single column
	#as well as the differences 
	st <- 1
	for(i in 1:noftrain){
		curseries <- trainingdata[i, !is.na(trainingdata[i, ])]
		curclass <- curseries[1]
		#standardize if necessary
		numseries <- as.numeric(curseries[2:length(curseries)])
		numseries <- (numseries-mean(numseries))/sd(numseries)
		en <- st+seriesLen[i]-2
		observations[st:en] <- numseries[2:length(numseries)]
		difference[st:en] <- diff(numseries)
		obsclass <- rep(curclass, seriesLen[i]-1)
		if(i == 1){
			allobsclass <- obsclass
		} else {
			allobsclass <- c(allobsclass, obsclass)
		}
		st <- en+1
	}
		
	timeindices <- unlist(lapply(seriesLen,function(x) c(2:x))) #create time indices
	#final train matrix stores class,time index,observation and consecutive difference
	finaltrain <- data.frame(Class = allobsclass, timeindices, observations, difference)
	ntrainobs <- seriesLen-1
	datatraintimeend <- proc.time()
	datatrainprepdur <- datatraintimeend-datatraintimestart
	datatrainprepdur <- datatrainprepdur[3]

	#algorithm for generating the distribution of symbols
	generatecodebook <- function(nodestatus, terminal, nofterminal, nofobservations) {
		if(!is.loaded("generate_codebook")) dyn.load("mts_functions64bit.dll")
		nofseries <- length(nofobservations)
		noftree <- ncol(terminal)
		nofnode <- nrow(nodestatus)
		total <- sum(nofobservations)
		nofentry <- nofseries*nofterminal*noftree
		out <- .C("generate_codebook", as.integer(as.matrix(nodestatus)), as.integer(nofnode), as.integer(noftree), as.integer(as.matrix(terminal)), as.integer(nofterminal), as.integer(nofobservations), as.integer(total), as.integer(nofseries), result = double(nofentry))
		return(out$result)
	}
	
	noftreelevels <- tuningParamLevels$noftreelevels
	nofnodelevels <- tuningParamLevels$nofnodelevels
	
	# source("parameter_selection_noparallel.r") 
	
	ntreeRFts <- 50
	t1 <- system.time({
		noftree <- 25
		OOB_error_rate_node <- array(0, length(nofnodelevels))
		for(nd in 1:length(nofnodelevels)){		# for each nnumber of trees (J_{ins})
			nofnode <- nofnodelevels[nd]
			RFins <- randomForest(as.matrix(finaltrain[, 2:ncol(finaltrain)]), factor(finaltrain[, 1]), ntree = noftree, maxnodes = nofnode)
			train_terminal <- attr(predict(RFins, finaltrain[, 2:ncol(finaltrain)], nodes = TRUE), "nodes")
			codetr <- matrix(generatecodebook(RFins$forest$nodestatus,train_terminal, nofnode, ntrainobs), noftrain, noftree*nofnode)	
			RFts <- randomForest(codetr, as.factor(classtrain), ntree = ntreeRFts)
			OOB_error_rate_node[nd] <- 1-sum(predict(RFts, type = "response") == as.factor(classtrain))/noftrain
		}
		#plot(nofnodelevels,OOB_error_rate_node,type="l",col="red",xlab="Alphabet size",ylab="Error rate",pch=4,lty=2)
		#points(nofnodelevels,OOB_error_rate_node,type="p",col="red",pch=4)
	})

	t2 <- system.time({
		OOB_error_rate <- array(0, length(noftreelevels))
		nofnode <- nofnodelevels[which.min(OOB_error_rate_node)]
		for(nd in 1:length(noftreelevels)){		# for each nnumber of trees (J_{ins})
			noftree <- noftreelevels[nd]
			RFins <- randomForest(finaltrain[, 2:ncol(finaltrain)], factor(finaltrain[, 1]), ntree = noftree, maxnodes = nofnode)
			train_terminal <- attr(predict(RFins, finaltrain[, 2:ncol(finaltrain)], nodes = TRUE), "nodes")
			codetr <- matrix(generatecodebook(RFins$forest$nodestatus, train_terminal, nofnode, ntrainobs), noftrain, noftree*nofnode)	
			RFts <- randomForest(codetr, as.factor(classtrain), ntree = ntreeRFts)
			OOB_error_rate[nd] <- 1-sum(predict(RFts, type = "response") == as.factor(classtrain))/noftrain
		}
		#plot(noftreelevels,OOB_error_rate,type="l",col="red",xlab="Number of trees",ylab="Error rate",pch=4,lty=2)
		#points(noftreelevels,OOB_error_rate,type="p",col="red",pch=4)
		noftree <- noftreelevels[which.min(OOB_error_rate)]
	})

	optParams = list(noftree = noftree, nofnode = nofnode)
	return(optParams)
	
}