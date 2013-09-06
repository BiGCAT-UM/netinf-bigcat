#!/usr/bin/Rscript

# HPN-DREAM Network Inference Challenge, 2013
# Gungor Budak (gungor.budak@maastrichtuniversity.nl)
# Department of Bioinformatics, Maastricht University

args <- commandArgs(TRUE)
midasFile <- args[1]
relations <- list(inhibitors = c("INH1_loLIG1", "INH1_hiLIG1",
				"INH1_loLIG2", "INH1_hiLIG2",
				"INH2_loLIG1", "INH2_hiLIG1",
				"INH2_loLIG2", "INH2_hiLIG2",
				"INH3_loLIG1", "INH3_hiLIG1",
				"INH3_loLIG2", "INH3_hiLIG2"),
		stimuli = 	c("loLIG1", "hiLIG1",
				"loLIG2", "hiLIG2",
				"loLIG1", "hiLIG1",
				"loLIG2", "hiLIG2",
				"loLIG1", "hiLIG1",
				"loLIG2", "hiLIG2"),
		targets = 	c("AB12_L1", "AB12_H1",
				"AB12_L2", "AB12_H2",
				"AB5_L1", "AB5_H1",
				"AB5_L2", "AB5_H2",
				"AB8_L1", "AB8_H1",
				"AB8_L2", "AB8_H2"),
		timePoints =	c(0, 1, 2, 4, 6, 10, 15, 30, 45, 60, 120))

midasReader <- function(midasFile) {

	midasData <- read.csv(file=midasFile, head=TRUE, sep=",")
	midasDataColNum <- length(midasData[1,])
	midasDataRowNum <- length(midasData[,1])

	dataList <- list(treatments = vector(), timePoints = vector())
	lastTreatment <- ""

	for (i in 1:midasDataRowNum) {

		treatment <- ""

		for (j in 1:midasDataColNum) {
		
			colName <- colnames(midasData[i,][j])
			colNameSplitted <- unlist(strsplit(colName, "[.]"))[2]

			if (grepl("TR", colName) == TRUE && grepl("CellLine", colName) == FALSE) {
				
				if (midasData[i,][[j]]) {

					if (treatment != "") { treatment <- paste(treatment, "_", sep="") }
					treatment <- paste(treatment, colNameSplitted, sep="")
					}
				} else if (grepl("DA", colName) == TRUE) {

					if (lastTreatment != treatment) {

						dataList[["timePoints"]] <- append(dataList[["timePoints"]], midasData[i,][[j]])
						}
					} else if (grepl("DV", colName) == TRUE) {

						if (lastTreatment != treatment) {

							average <- (midasData[i,][[j]] + midasData[i+1,][[j]] + midasData[i+2,][[j]])/3
							dataList[[colNameSplitted]] <- append(dataList[[colNameSplitted]], average)
							}
						}
			} # End of each column

			if (lastTreatment != treatment) {
				dataList[["treatments"]] <- append(dataList[["treatments"]], treatment)
				lastTreatment <- treatment
				}
		} # End of each row

	return(dataList)
	}

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

					if (is.na(dataList[[j]][[k]]) == FALSE) {
						optimizedData[[intervention]][[names(dataList[j])]][[counter1]] <- dataList[[j]][[k]]
						counter1 <- counter1 + 1
						}
					}

				if (dataList$treatments[[k]] == noIntervention) {

					if (is.na(dataList[[j]][[k]]) == FALSE) {
						optimizedData[[noIntervention]][[names(dataList[j])]][[counter2]] <- dataList[[j]][[k]]
						counter2 <- counter2 + 1
						}
					}	
				}
			}
		}
	
	return(optimizedData)
	}

edgeScorer <- function(optimizedData) {
	
	maxDifference <- 0
	scores <- list()

	# Find the maximum value of differences to normalize scores

	for (i in 1:length(relations$targets)) {

		for (j in 1:length(optimizedData[[relations$inhibitors[[i]]]])) {
			
			currentDifference <- abs(optimizedData[[relations$inhibitors[[i]]]][[j]] - optimizedData[[relations$stimuli[[i]]]][[j]])
			if (currentDifference > maxDifference) {

				maxDifference <- currentDifference
				}
			}
		}

	# Normalize scores

	for (i in 1:length(relations$targets)) {

		for (j in 1:length(optimizedData[[relations$inhibitors[[i]]]])) {

			currentDifference <- abs(optimizedData[[relations$inhibitors[[i]]]][[j]] - optimizedData[[relations$stimuli[[i]]]][[j]])
			ratio <- currentDifference/maxDifference
			target <- relations$targets[[i]]
			antibody <- names(optimizedData[[relations$inhibitors[[i]]]][j])
			scoreId <- paste(target, antibody, sep="_")
			scores[[scoreId]] <- ratio
			}
		}

	return(scores)
	}

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

		for (j in 1:length(optimizedData[[relations$inhibitors[[i]]]])) {

			comparison <- (optimizedData[[relations$inhibitors[[i]]]][[j]] > optimizedData[[relations$stimuli[[i]]]][[j]])
			relation <- if (comparison) "-1" else "1"
			target <- relations$targets[[i]]
			antibody <- names(optimizedData[[relations$inhibitors[[i]]]][j])
			scoreId <- paste(target, antibody, sep="_")

			targetOriginal <- unlist(strsplit(target, "[_]"))[1]
			scoreIdOriginal <- paste(targetOriginal, antibody, sep="_")

			if (scoresOptimized[[scoreId]] > 0.1) {

				# Check if the edge is already in results
				theMatch <- match(scoreIdOriginal, results$edge)
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
						results$edge <- append(results$edge, paste(targetOriginal, antibody, sep="_"))
						results$relation <- append(results$relation, relation)
						results$score <- append(results$score, scoresOptimized[[scoreId]])
						}
				}
			}
		}

	# Output results

	sink("output/BiGCaT-Network-Insilico.eda")
	for (i in 1:length(results$edge)) {

		cat(paste(
			unlist(strsplit(results$edge[[i]], "[_]"))[1],
			"\t(",
			results$relation[[i]],
			")\t",
			unlist(strsplit(results$edge[[i]], "[_]"))[2],
			"\t=\t",
			results$score[[i]],
			"\n",
			sep=""))
		}
	sink()

	sink("output/BiGCaT-Network-Insilico.sif")
	for (i in 1:length(results$edge)) {

		cat(paste(
			unlist(strsplit(results$edge[[i]], "[_]"))[1],
			"\t",
			results$relation[[i]],
			"\t",
			unlist(strsplit(results$edge[[i]], "[_]"))[2],
			"\n",
			sep=""))
		}
	sink()
	}

dataList <- midasReader(midasFile)
optimizedData <- dataOptimizer(dataList)
networkInferrer(optimizedData)
