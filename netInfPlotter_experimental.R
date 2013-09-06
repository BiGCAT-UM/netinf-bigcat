#!/usr/bin/Rscript

# HPN-DREAM Network Inference Challenge, 2013
# Gungor Budak (gungor.budak@maastrichtuniversity.nl)
# Department of Bioinformatics, Maastricht University

args <- commandArgs(TRUE)
cellLine <- args[1]
childAntibody <- args[2]
noIntervention <- args[3]
intervention <- args[4]
midasFilePath <- paste("MD_", cellLine, "_main.csv", sep="")

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

# This function removes duplicates by
# taking their average and assign as one

duplicateHandler <- function(experiments) {

	lastTreatment <- ""
	numberOfRows <- nrow(experiments)
	dataList <- list(
		treatments = vector(),
		timePoints = vector(),
		dataValues = vector(),
		orders = vector()
		)
	count <- 1

	for (i in 1:numberOfRows) {
		
		if (lastTreatment != experiments$treatments[[i]]) {

			dataList$treatments[[count]] <- toString(experiments$treatments[[i]])
			dataList$timePoints[[count]] <- experiments$timePoints[[i]]
			dataList$dataValues[[count]] <- experiments$dataValues[[i]]
			dataList$orders[[count]] <- experiments$orders[[i]]

			lastTreatment <- experiments$treatments[[i]]
			count <- count+1
			} else {

				dataList$dataValues[[which(dataList$orders == experiments$orders[[i]])]] <- (experiments$dataValues[[i]] + dataList$dataValues[[which(dataList$orders == experiments$orders[[i]])]])/2
				}
		}

	clearedExperiments <- data.frame(treatments=dataList$treatments, timePoints=dataList$timePoints, dataValues=dataList$dataValues, orders=dataList$orders)

	return(clearedExperiments)
	}

# This function knows what to have in the data
# and checks data according to that
# if there is missing data, it put NA for data value

fashionChecker <- function(experiments) {

	numberOfRows <- nrow(experiments)
	dataList <- list(
		treatments = vector(),
		timePoints = vector(),
		dataValues = vector()
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

		if (treatment != experiments$treatments[[l]]) {

			dataList$treatments[[z]] <- treatment
			dataList$timePoints[[z]] <- experiments$timePoints[[l]]
			dataList$dataValues[[z]] <- NA
			} else {

				dataList$treatments[[z]] <- toString(experiments$treatments[[l]])
				dataList$timePoints[[z]] <- experiments$timePoints[[l]]
				dataList$dataValues[[z]] <- experiments$dataValues[[l]]
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

	# completeExperiments <- data.frame(treatments=dataList$treatments, timePoints=dataList$timePoints, dataValues=dataList$dataValues)

	return(dataList)
	}

# Function for estimation missing data points
# only when they are single and in the middle

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

# This is a function that read MIDAS file into
# a data frame and optimize the data

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

						childAntibodyOriginal <- paste("DV.", childAntibody, sep="")
						if (colName == childAntibodyOriginal) {
							dataList[[childAntibody]] <- append(dataList[[childAntibody]], midasData[i,][[j]])
							}
						}
			} # End of each column

			dataList$treatments <- append(dataList$treatments, treatment)

		} # End of each row

	
	
	experiments <- data.frame(treatments=dataList$treatments, timePoints=dataList$timePoints, dataValues=dataList[[childAntibody]])

	# Get rid of 0 time point
	experiments <- subset(experiments, timePoints!=0)

	experiments <- orderData(experiments)
	clearedExperiments <- duplicateHandler(experiments)
	completeExperiments <- fashionChecker(clearedExperiments)
	estimatedExperiments <- estimateDataPoints(completeExperiments)

	estimatedExperiments <- do.call(cbind.data.frame, estimatedExperiments)

	return(estimatedExperiments)
	} # End of midasReader function

# Plot the data for child node and target node

dataPlotter <- function(dataList, noIntervention, intervention) {

	numberOfRows <- nrow(dataList)
	xAxis <- vector()
	yAxisNoIntervention <- vector()
	yAxisIntervention <- vector()

	for (i in 1:numberOfRows) {
		
		# Get data for when there is NO intervention on the target
		if (dataList$treatments[[i]] == noIntervention) {

			# Create x-axis vector
			xAxis <- append(xAxis, dataList$timePoints[[i]])
			yAxisNoIntervention <- append(yAxisNoIntervention, dataList$dataValues[[i]])

			# Get data for when there is intervention on the target
			} else if (dataList$treatments[[i]] == intervention) {
				
				yAxisIntervention <- append(yAxisIntervention, dataList$dataValues[[i]])
				}
		}

	# Do not plot the graph with missing data
	if ((NA %in% yAxisNoIntervention) || (NA %in% yAxisIntervention)) {
		print("Data value missing cannot plot the graph")
		} else {
			yAxisMax <- range(0, yAxisNoIntervention, yAxisIntervention)
			plot(xAxis, yAxisNoIntervention, xlab="Time Points", ylab=paste(gsub("DV.", "", childAntibody), "Data", sep=" "),
				type = "b", col = "blue", ylim = yAxisMax)
			lines(xAxis, yAxisIntervention, type = "b", pch = 22, col="red", lty = 2)
			legend("bottomright", c("No Intervention", "Intervention"),
				cex=0.8, col = c("blue", "red"), pch=21:22, lty=1:2, lwd=2, bty="n");
			title(main=paste(gsub("DV.", "", childAntibody), " Expression Data ", "(Cell line: ", cellLine, " & ", "Stimulus: ", noIntervention, ")", sep=""),
				col.main="red", font.main=4)
			}
	} # End of dataPlotter function

dataList <- midasReader(midasFilePath, childAntibody)
dataPlotter(dataList, noIntervention, intervention)
