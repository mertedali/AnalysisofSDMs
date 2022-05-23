predict.SMTS <- function(newdata, modelPath){

	# setwd("C:/Users/Mert/Desktop/R Package/R Package v2")
	
	# source("newdata_preparation.r")
	
	nofnew <- nrow(newdata)
	seriesLen <- apply(newdata, 1, function(x) sum(!is.na(x)))
	observations <- array(0, sum(seriesLen)-nofnew)
	difference <- array(0, sum(seriesLen)-nofnew)
	st <- 1
	for(i in 1:nofnew){
		curseries = newdata[i, !is.na(newdata[i, ])]
		#standardize if necessary
		numseries <- as.numeric(curseries)
		numseries <- (numseries-mean(numseries))/sd(numseries)
		en <- st+seriesLen[i]-2
		observations[st:en] <- numseries[2:length(numseries)]
		difference[st:en] <- diff(numseries)
		st <- en+1
	}
		
	timeindices <- unlist(lapply(seriesLen, function(x) c(2:x)))	#create time indices
	finalnew <- data.frame(timeindices, observations, difference)
	nnewobs <- seriesLen-1

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
	
	model <- readRDS(modelPath)

	new_terminal <- attr(predict(model$RFins, as.matrix(finalnew), nodes = TRUE), "nodes")
	codenew <- matrix(generatecodebook(model$RFins$forest$nodestatus, new_terminal, model$nofnode, nnewobs), nofnew, model$noftree*model$nofnode)
	
	predicted_prob <- predict(model$RFts, codenew, type = "prob")
	colnames(predicted_prob) <- model$classInfo$x
	predicted_class <- predict(model$RFts, codenew, type = "class")
	
	predicted_class_new <- array(NA, length(predicted_class))
	for(j in 1:nrow(model$classInfo)){
		predicted_class_new[which(predicted_class == model$classInfo$ID[j])] <- as.character(model$classInfo$x[j])
	}
	 
	return(list(classPred = as.factor(predicted_class_new), probVals = predicted_prob))

}