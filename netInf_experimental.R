#!/usr/bin/Rscript

# HPN-DREAM Network Inference Challenge, 2013
# Gungor Budak (gungor.budak@maastrichtuniversity.nl)
# Department of Bioinformatics, Maastricht University

# netInf_experimental.R - Network Inference Using Expermental Data from the Challenge

# This script works R 2.14 and above and requires simp.R from StreamMetabolism package in the same directory
# takes cell line and stimulus as arguments; cell line is used to locate necessary data file and stimulus
# is used to characterize the data according to intervention and no-intervention studies.

args <- commandArgs(TRUE)
cellLine <- args[1]
stimulus <- args[2]
midasFilePath <- paste("MD_", cellLine, "_main.csv", sep="")

# relations data list holds relations between stimuli / inhibitors and their target proteins
# Note that inhibitor - target protein relation information is obtained from literature

relations <- list(
		inhibitors=c(paste(stimulus, "GSK690693", sep="__"), paste(stimulus, "GSK690693_GSK1120212", sep="__"), paste(stimulus, "PD173074", sep="__")),
		stimuli = c(stimulus, stimulus, stimulus),
		targets = c("AKT", "AKT__MEK", "FGFR1__FGFR3"),
		timePoints = c(5, 15, 30, 60, 120, 240))

# orderData function rearranges the characterized data points according to a fashion so that
# later we can see if there is any missing data point and which treatment is missing

orderData <- function(experiments) {

	numberOfRows <- nrow(experiments)
	orders <- vector()
	dataFashion <- c("FGF1__GSK690693",
			"Insulin__GSK690693",
			"FGF1__GSK690693_GSK1120212",
			"Insulin__GSK690693_GSK1120212",
			"FGF1__PD173074",
			"Insulin__PD173074",
			"FGF1",
			"Insulin",
			"EGF__GSK690693",
			"IGF1__GSK690693",
			"EGF__GSK690693_GSK1120212",
			"IGF1__GSK690693_GSK1120212",
			"EGF__PD173074",
			"IGF1__PD173074",
			"EGF",
			"IGF1",
			"HGF__GSK690693",
			"Serum__GSK690693",
			"HGF__GSK690693_GSK1120212",
			"Serum__GSK690693_GSK1120212",
			"HGF__PD173074",
			"Serum__PD173074",
			"HGF",
			"Serum",
			"NRG1__GSK690693",
			"PBS__GSK690693",
			"NRG1__GSK690693_GSK1120212",
			"PBS__GSK690693_GSK1120212",
			"NRG1__PD173074",
			"PBS__PD173074",
			"NRG1",
			"PBS")

	order <- 1
	tPConstant <- 0
	tP <- 5

	for (i in 1:numberOfRows) {

		if (experiments$timePoints[[i]] != tP) {
			tPConstant <- tPConstant+1
			tP <- experiments$timePoints[[i]]
			}

		order <- (which(dataFashion==experiments$treatment[[i]])) + (32*tPConstant)
		orders <- append(orders, order)
		}

	experiments <- transform(experiments, orders = orders)

	experiments <- experiments[order(experiments$orders),]

	return(experiments)
	}

# duplicateHandler looks for replicates and if there are any, it averages them into a single data point

duplicateHandler <- function(experiments) {

	numberOfRows <- nrow(experiments)
	numberOfCols <- ncol(experiments)
	dataList <- list(
		orders = vector(),
		treatments = vector(),
		timePoints = vector()
		)

	for (i in 1:numberOfRows) {

		for (j in 1:numberOfCols) {
		
			colName <- toString(names(experiments[j]))
			pos <- which(dataList$orders == experiments$orders[[i]])

			if (length(pos)) {
				
				if (colName != "orders" && colName != "treatments" && colName != "timePoints") {
					dataList[[colName]][[pos]] <- (dataList[[colName]][[pos]] + experiments[[colName]][[i]])/2
					}
				} else {

					if (colName == "orders") {
						dataList$orders <- append(dataList$orders, experiments[i,][[j]])
						} else if (colName == "treatments") {
							dataList$treatments <- append(dataList$treatments, toString(experiments[i,][[j]]))
							} else if (colName == "timePoints") {
								dataList$timePoints <- append(dataList$timePoints, experiments[i,][[j]])
								} else {
									dataList[[colName]] <- append(dataList[[colName]], experiments[i,][[j]])
									}
					}


			}
		}

	clearedExperiments <- do.call(cbind.data.frame, dataList)

	return(clearedExperiments)
	}

# fashinChecker knows what to have in the data and checks data according to that
# if there is missing data. If so, it puts NA for the missing data value

fashionChecker <- function(experiments) {

	numberOfRows <- nrow(experiments)
	numberOfCols <- ncol(experiments)
	dataList <- list(
		treatments = vector(),
		timePoints = vector()
		)

	dataFashion <- list(
		inhibitors = c("GSK690693", "GSK690693_GSK1120212", "PD173074", ""),
		stimuli = c("FGF1", "Insulin", "EGF", "IGF1", "HGF", "Serum", "NRG1", "PBS")
		)

	j <- 1
	k <- 1
	l <- 1
	i <- 1
	z <- 1

	while (l <= numberOfRows) {
		
		if (i == 33) {
			i <- 1
			}

		check = ((i == 3) || (i == 5) || (i == 7) ||
			(i == 11) || (i == 13) || (i == 15) ||
			(i == 19) || (i == 21) || (i == 23) ||
			(i == 27) || (i == 29) || (i == 31))

		if (check) {
			k <- k+1
			}

		if ((i %% 2) != 0) {
		
			if (dataFashion$inhibitors[[k]] != "") {	
				treatment <- paste(dataFashion$stimuli[[j]], dataFashion$inhibitors[[k]], sep="__")
				} else {
					treatment <- dataFashion$stimuli[[j]]
					}
			} else {

				if (dataFashion$inhibitors[[k]] != "") {	
					treatment <- paste(dataFashion$stimuli[[j+1]], dataFashion$inhibitors[[k]], sep="__")
					} else {
						treatment <- dataFashion$stimuli[[j+1]]
						}
				}
		
		y <- 1

		while (y <= numberOfCols) {

			colName <- toString(names(experiments[y]))

			if (treatment != experiments$treatments[[l]]) {

				dataList$treatments[[z]] <- treatment
				dataList$timePoints[[z]] <- experiments$timePoints[[l]]
				if (colName != "orders" && colName != "treatments" && colName != "timePoints") {
					dataList[[colName]][[z]] <- NA
					}
				inc <- 0
				} else {

					dataList$treatments[[z]] <- toString(experiments$treatments[[l]])
					dataList$timePoints[[z]] <- experiments$timePoints[[l]]
					if (colName != "orders" && colName != "treatments" && colName != "timePoints") {
						dataList[[colName]][[z]] <- experiments[[colName]][[l]]
						}
					inc <- 1
					}
			y <- y+1
			}

		if (inc) {
			l <- l+1
			}

		if ((i %% 8) == 0) {
			j <- j+2
			k <- 1
			}

		if ((i %% 32) == 0) {
			j <- 1
			}

		i <- i+1
		z <- z+1

		}

	return(dataList)
	}


# estimateDataPoints looks for NAs, and if they are not at the beginning or at the end
# and if there are not more than one NA in a row, it estimates that data point
# by taking the average of neighboring data points in time series data

estimateDataPoints <- function(experiments) {

	numberOfRows <- length(experiments[[1]])
	numberOfCols <- length(experiments)

	for (i in 1:numberOfRows) {
		for (j in 1:numberOfCols) {

			if (experiments$timePoints[[i]] != 5 && experiments$timePoints[[i]] != 240) {

				if (is.na(experiments[[j]][[i]]) == TRUE) {

					if (is.na(experiments[[j]][[i-32]]) == FALSE && is.na(experiments[[j]][[i+32]]) == FALSE) {
						experiments[[j]][[i]] <- ( experiments[[j]][[i-32]] + experiments[[j]][[i+32]] ) / 2
						}

					}
				}
			}
		}

	return(experiments)
	}

# midasReader is a function that reads MIDAS file into
# a data frame and prepares the data for further analyses by rearranging
# handling replicates and estimating missing data points

midasReader <- function(midasFilePath, childAntibody) {

	# Read data into data object
	midasData <- read.csv(file=midasFilePath, head=TRUE, sep=",")
	midasDataColNum <- length(midasData[1,])
	midasDataRowNum <- length(midasData[,1])

	# List to hold data
	dataList <- list(
		treatments = vector(),
		timePoints = vector()
		)

	# This loop goes over each row in MIDAS file
	for (i in 1:midasDataRowNum) {

		treatment <- ""

		for (j in 1:midasDataColNum) {

			colName <- colnames(midasData[i,][j])
	
			if (grepl("TR.", colName) == TRUE && grepl("CellLine", colName) == FALSE) {

				if (midasData[i,][[j]]) {

					if (treatment != "") {

						treatment <- paste(treatment, "__", sep="")
						}
					colNameSplitted <- unlist(strsplit(colName, "[.]"))[2]
					treatment <- paste(treatment, colNameSplitted, sep="")
					}

				} else if (grepl("DA.", colName) == TRUE) {
					
					dataList$timePoints <- append(dataList$timePoints, midasData[i,][[j]])
					} else if (grepl("DV.", colName) == TRUE) {

							colNameSplitted <- gsub("DV.", "", colName)
							colNameSplitted <- gsub("[.]", "", colNameSplitted)
							dataList[[colNameSplitted]] <- append(dataList[[colNameSplitted]], midasData[i,][[j]])
						}
			} # End of each column

			dataList$treatments <- append(dataList$treatments, treatment)

		} # End of each row

	
	
	experiments <- do.call(cbind.data.frame, dataList)

	# Get rid of 0 time point
	experiments <- subset(experiments, timePoints!=0)

	experiments <- orderData(experiments)
	clearedExperiments <- duplicateHandler(experiments)
	completeExperiments <- fashionChecker(clearedExperiments)
	estimatedExperiments <- estimateDataPoints(completeExperiments)

	return(estimatedExperiments)
	} # End of midasReader function

# dataOptimizer characterizes that data according to intervention
# and no intervention studies

dataOptimizer <- function(dataList) {

	optimizedData <- list()

	for (i in 1:length(relations$targets)) {

		intervention <- relations$inhibitors[i]
		noIntervention <- relations$stimuli[i]
		optimizedData[[intervention]] <- list()
		optimizedData[[noIntervention]] <- list()

		# First two vectors have treatments and time points
		for (j in 3:length(dataList)) {
			# For each node

			counter1 <- 1
			counter2 <- 1
			for (k in 1:length(dataList$treatments)) {
				# For each time point
				
				if (dataList$treatments[[k]] == intervention) {

					optimizedData[[intervention]][[names(dataList[j])]][[counter1]] <- dataList[[j]][[k]]
					counter1 <- counter1 + 1
					}

				if (dataList$treatments[[k]] == noIntervention) {

					optimizedData[[noIntervention]][[names(dataList[j])]][[counter2]] <- dataList[[j]][[k]]
					counter2 <- counter2 + 1
					}		
				}
			}
		}

	return(optimizedData)
	}

# edgeScorer scores the edges by finding max difference AUC
# and divides all to max to have scores between 0 and 1
# this shows have diverged intervention and no-intervention from each other
# in terms of protein expression

edgeScorer <- function(optimizedData) {
	
	maxDifference <- 0
	scores <- list()

	# Find the maximum value of differences to normalize scores

	for (i in 1:length(relations$targets)) {

		if (any(is.na(optimizedData[[relations$inhibitors[[i]]]])) == FALSE) {

			for (j in 1:length(optimizedData[[relations$inhibitors[[i]]]])) {
			
				currentDifference <- abs(optimizedData[[relations$inhibitors[[i]]]][[j]] - optimizedData[[relations$stimuli[[i]]]][[j]])
				if (currentDifference > maxDifference) {

					maxDifference <- currentDifference
					}
				}
			}
		}

	# Normalize scores

	for (i in 1:length(relations$targets)) {

		if (any(is.na(optimizedData[[relations$inhibitors[[i]]]])) == FALSE) {

			for (j in 1:length(optimizedData[[relations$inhibitors[[i]]]])) {

				currentDifference <- abs(optimizedData[[relations$inhibitors[[i]]]][[j]] - optimizedData[[relations$stimuli[[i]]]][[j]])
				ratio <- currentDifference/maxDifference
				target <- relations$targets[[i]]
				antibody <- names(optimizedData[[relations$inhibitors[[i]]]][j])
				scoreId <- paste(target, antibody, sep="___")
				scores[[scoreId]] <- ratio
				}
			}
		}

	# Separate edges with two targets

	finalScores <- list()

	for (i in 1:length(scores)) {

		scoreName <- names(scores[i])
		if (grepl("[a-zA-Z0-9]+(_){2}[a-zA-Z0-9]+", scoreName, perl=TRUE) == TRUE) {

			targets <- unlist(strsplit(scoreName, "___", fixed=TRUE))[1]
			antibody <- unlist(strsplit(scoreName, "___", fixed=TRUE))[2]

			target1 <- unlist(strsplit(targets, "__", fixed=TRUE))[1]
			target2 <- unlist(strsplit(targets, "__", fixed=TRUE))[2]

			newScoreId1 <- paste(target1, antibody, sep="___")
			newScoreId2 <- paste(target2, antibody, sep="___")

			finalScores[[newScoreId1]] <- scores[[i]]
			finalScores[[newScoreId2]] <- scores[[i]]
			} else {

				if (length(finalScores[[scoreName]])) {
					
					if (scores[[i]] > finalScores[[scoreName]]) {
						finalScores[[scoreName]] <- scores[[i]]
						}
					} else {
						finalScores[[scoreName]] <- scores[[i]]
						}
				
				}
		}

	return(scores)
	}
	
# networkInferrer gets the final form of data and infers the network
# and produces SIF and EDA files of the inferred network

networkInferrer <- function(optimizedData) {

	source("simp.R")
	results <- list(edge=vector(), relation=vector(), score=vector())
	

	# Integrate time points for each condition

	for (i in 1:length(optimizedData)) {

		for (j in 1:length(optimizedData[[i]])) {

			optimizedData[[i]][[j]] <- simp(x=relations$timePoints, y=optimizedData[[i]][[j]])
			}
		}

	# Get scores of reoptimized (by integrating) data

	scoresOptimized <- edgeScorer(optimizedData)

	# Infer network

	for (i in 1:length(relations$targets)) {

		# We need to make sure that there is no missing data
		if (any(is.na(optimizedData[[relations$inhibitors[[i]]]])) == FALSE) {

			for (j in 1:length(optimizedData[[relations$inhibitors[[i]]]])) {

				comparison <- (optimizedData[[relations$inhibitors[[i]]]][[j]] > optimizedData[[relations$stimuli[[i]]]][[j]])
				relation <- if (comparison) "-1" else "1"
				target <- relations$targets[[i]]
				antibody <- names(optimizedData[[relations$inhibitors[[i]]]][j])
				scoreId <- paste(target, antibody, sep="___")

				# Get edges scored larger than 0.1
				if (scoresOptimized[[scoreId]] > 0.1) {

					# Check if the edge is already in results
					theMatch <- match(scoreId, results$edge)
					if (is.na(theMatch) == FALSE) {

						# Check if the edge present has low score
						if (scoresOptimized[[scoreId]] > results$score[[theMatch]]) {

							# Modification of previously added edge
							# in terms of relation and score based on their scores
							results$relation[[theMatch]] <- relation
							results$score[[theMatch]] <- scoresOptimized[[scoreId]]
							}
						} else {

							# Addition of new edge to results
							results$edge <- append(results$edge, paste(target, antibody, sep="___"))
							results$relation <- append(results$relation, relation)
							results$score <- append(results$score, scoresOptimized[[scoreId]])
							}
					}
				}
			}
		}

	# Separate edges with two targets

	finalResults <- list(edge=vector(), relation=vector(), score=vector())

	for (i in 1:length(results$edge)) {

		resultName <- results$edge[[i]]
		
		if (grepl("[a-zA-Z0-9]+(_){2}[a-zA-Z0-9]+", resultName, perl=TRUE) == TRUE) {

			targets <- unlist(strsplit(resultName, "___", fixed=TRUE))[1]
			antibody <- unlist(strsplit(resultName, "___", fixed=TRUE))[2]

			target1 <- unlist(strsplit(targets, "__", fixed=TRUE))[1]
			target2 <- unlist(strsplit(targets, "__", fixed=TRUE))[2]

			newResultName1 <- paste(target1, antibody, sep="___")
			newResultName2 <- paste(target2, antibody, sep="___")

			pos <- which(finalResults$edge == newResultName1)
			if (length(pos)) {

				if (results$score[[i]] > finalResults$score[[pos]]) {
					finalResults$score[[pos]] <- results$score[[i]]
					finalResults$relation[[pos]] <- results$relation[[i]]
					}
				} else {
					finalResults$edge <- append(finalResults$edge, newResultName1)
					finalResults$relation <- append(finalResults$relation, results$relation[[i]])
					finalResults$score <- append(finalResults$score, results$score[[i]])
					}

			pos <- which(finalResults$edge == newResultName2)
			if (length(pos)) {

				if (results$score[[i]] > finalResults$score[[pos]]) {
					finalResults$score[[pos]] <- results$score[[i]]
					finalResults$relation[[pos]] <- results$relation[[i]]
					}
				} else {
					finalResults$edge <- append(finalResults$edge, newResultName2)
					finalResults$relation <- append(finalResults$relation, results$relation[[i]])
					finalResults$score <- append(finalResults$score, results$score[[i]])
					}
			} else {

				finalResults$edge <- append(finalResults$edge, resultName)
				finalResults$relation <- append(finalResults$relation, results$relation[[i]])
				finalResults$score <- append(finalResults$score, results$score[[i]])
				}
		}

	# Name targets

	outputResults <- list(edge=vector(), relation=vector(), score=vector())

	for (i in 1:length(finalResults$edge)) {

		resultName <- finalResults$edge[[i]]
		target <- unlist(strsplit(resultName, "___", fixed=TRUE))[1]
		antibody <- unlist(strsplit(resultName, "___", fixed=TRUE))[2]

		if (target == "AKT") {

			newResultName1 <- paste("AKT_pS473", antibody, sep="___")
			outputResults$edge <- append(outputResults$edge, newResultName1)
			outputResults$relation <- append(outputResults$relation, finalResults$relation[[i]])
			outputResults$score <- append(outputResults$score, finalResults$score[[i]])

			newResultName2 <- paste("AKT_pT308", antibody, sep="___")
			outputResults$edge <- append(outputResults$edge, newResultName2)
			outputResults$relation <- append(outputResults$relation, finalResults$relation[[i]])
			outputResults$score <- append(outputResults$score, finalResults$score[[i]])
			} else if (target == "MEK") {

				newResultName1 <- paste("MAPK_pT202_Y204", antibody, sep="___")
				outputResults$edge <- append(outputResults$edge, newResultName1)
				outputResults$relation <- append(outputResults$relation, finalResults$relation[[i]])
				outputResults$score <- append(outputResults$score, finalResults$score[[i]])

				newResultName2 <- paste("MEK1_pS217_S221", antibody, sep="___")
				outputResults$edge <- append(outputResults$edge, newResultName2)
				outputResults$relation <- append(outputResults$relation, finalResults$relation[[i]])
				outputResults$score <- append(outputResults$score, finalResults$score[[i]])
				} else {

					outputResults$edge <- append(outputResults$edge, resultName)
					outputResults$relation <- append(outputResults$relation, finalResults$relation[[i]])
					outputResults$score <- append(outputResults$score, finalResults$score[[i]])
					}

		}

	# Output results

	# Output edge scores
	sink(paste("output/BiGCaT-", cellLine, "-", stimulus, "-Network.eda", sep=""))
	for (i in 1:length(outputResults$edge)) {

		# We're asked to exclude these phopshoproteins
		if (grepl("TAZ_pS89", outputResults$edge[[i]]) == FALSE && grepl("FOXO3a_pS318_S321", outputResults$edge[[i]]) == FALSE) {

			cat(paste(
				unlist(strsplit(outputResults$edge[[i]], "___", fixed=TRUE))[1],
				"\t(",
				outputResults$relation[[i]],
				")\t",
				unlist(strsplit(outputResults$edge[[i]], "___", fixed=TRUE))[2],
				"\t=\t",
				outputResults$score[[i]],
				"\n",
				sep=""))
			}
		}
	sink()

	# Output network
	sink(paste("output/BiGCaT-", cellLine, "-", stimulus, "-Network.sif", sep=""))

	for (i in 1:length(outputResults$edge)) {

		# We're asked to exclude these phopshoproteins
		if (grepl("TAZ_pS89", outputResults$edge[[i]]) == FALSE && grepl("FOXO3a_pS318_S321", outputResults$edge[[i]]) == FALSE) {

			cat(paste(
				unlist(strsplit(outputResults$edge[[i]], "___", fixed=TRUE))[1],
				"\t",
				outputResults$relation[[i]],
				"\t",
				unlist(strsplit(outputResults$edge[[i]], "___", fixed=TRUE))[2],
				"\n",
				sep=""))
			}
		}
	sink()
	}

dataList <- midasReader(midasFilePath)
optimizedData <- dataOptimizer(dataList)
networkInferrer(optimizedData)
