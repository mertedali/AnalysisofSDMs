train.SMTS <- function(trainingdata, classes, tuningParams = list(noftree = 50, nofnode = 10), params = list(maxiter = 20, noftree_step = 50, tolerance = 0.05), saveModel = TRUE, savePath = getwd()){

	# noftreelevels: number of trees for symbol generation J_{ins}
	# nofnodelevels: alphabet size R
	# maxiter: the maximum number of iterations for training trees
	# tolerance: the number of trees to train per iteration
	
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
	
	noftree <- tuningParams$noftree
	nofnode <- tuningParams$nofnode
	
	RFins <- randomForest(as.matrix(finaltrain[, 2:ncol(finaltrain)]), factor(finaltrain[, 1]), ntree = noftree, maxnodes = nofnode)

	#get terminal node ids	
	train_terminal <- attr(predict(RFins, as.matrix(finaltrain[, 2:ncol(finaltrain)]), nodes = TRUE), "nodes")	
	codetr <- matrix(generatecodebook(RFins$forest$nodestatus, train_terminal, nofnode, ntrainobs), noftrain, noftree*nofnode)

	noftree_step <- params$noftree_step
	tolerance <- params$tolerance
	maxiter <- params$maxiter
	
	RFts <- randomForest(codetr, factor(classtrain), ntree = noftree_step)
	prev_OOBerror <- 1; cur_OOBerror <- 1-sum(predict(RFts, type = "response") == classtrain)/noftrain
	iter <- 1
	while(iter < 20 && cur_OOBerror < (1-tolerance)*prev_OOBerror){    
		prev_OOBerror <- cur_OOBerror
		RFtsmid <- randomForest(codetr, factor(classtrain), ntree = noftree_step, nodesize = 2)
		RFts <- combine(RFts, RFtsmid)
		cur_OOBerror <- 1-sum(predict(RFts, type = "response") == classtrain)/noftrain
		iter <- iter+1
	}	   
	OOB_error <- 1-sum(predict(RFts, type = "response") == classtrain)/noftrain
	
	model <- list(RFins = RFins, RFts = RFts, noftree = noftree, nofnode = nofnode, classInfo = uniqueclass)
	
	if(saveModel == TRUE){
		saveRDS(model, paste0(savePath, "/model.rds"))
		sprintf("Your model is saved to the following path: %s", savePath)
		sprintf("OOB error is %f", OOB_error)
	} else {
		return(model)
		sprintf("OOB error is %f", OOB_error)
	}
}