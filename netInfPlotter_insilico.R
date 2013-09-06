#!/usr/bin/Rscript

# HPN-DREAM Network Inference Challenge, 2013
# Gungor Budak (gungor.budak@maastrichtuniversity.nl)
# Department of Bioinformatics, Maastricht University

args <- commandArgs(TRUE)
midasFilePath <- args[1]
childAntibody <- args[2]
noIntervention <- args[3]
intervention <- args[4]


# This is a function that read MIDAS file into
# a data frame and optimize the data
midasReader <- function(midasFilePath, childAntibody){

	# Read data into data object
	proteomicsData <- read.csv(file=midasFilePath, head=TRUE, sep=",")

	# Vectors to hold data
	treatments <- vector()
	timePoints <- vector()
	data <- vector()
	
	#Because antibodies are given with "DV." in front default in MIDAS
	childAntibody <- paste("DV.", childAntibody, sep="")

	lastlyAddedTreatment <- ""

	# This is added specific to dataset, needs to be done dynamically
	numberOfTreatments <- 7 + 1
	numberOfExperiments <- length(proteomicsData[,1])

	# This loop goes over each row in MIDAS file
	for (i in 1:numberOfExperiments) {

		treatment <- ""

		# Create the vector of treatments (inhibitors/stimuli)
		# Starts from 2 because the first column is for cell line
		# This loop goes over each column of treatments (except cell line) in MIDAS file
		for (j in 2:numberOfTreatments) {
	
			# If the treatment is present ... (denoted as a "1" in MIDAS file)
			if (proteomicsData[i,][[j]]) {

				# Remove "TR." from the names
				permenantName <- unlist(strsplit(colnames(proteomicsData[i,][j]), "[.]"))[2]
				# If another treatment is present concatanate them with a plus sign
				if (treatment != "") { treatment <- paste(treatment, "+", sep="") }
				# There is extra "i" at the end of inhibitors remove it, too
				if (grepl("INH", permenantName)) { permenantName <- gsub("i", "", permenantName) }
				
				treatment <- paste(treatment, permenantName, sep="")
				}
			} # End of vertical loop for getting treatments

		# Don't add the same treatment twice
		if (lastlyAddedTreatment != treatment) {
			# This is the first row from the top for each treatment
			treatments <- append(treatments, treatment)
			lastlyAddedTreatment <- treatment
			# It's okay to take time point from the first row
			timePoints <- append(timePoints, proteomicsData$DA.ALL[[i]])

			# Since we're in the first row, the rest is +1 and +2
			# Taking average of three duplicates for the experiment
			data <- append(data, ((	proteomicsData[[childAntibody]][[i]] + 
						proteomicsData[[childAntibody]][[i+1]] +
						proteomicsData[[childAntibody]][[i+2]] ) / 3 ))
			}
		} # End of the horizontal loop for getting each experiment

	# Final data frame to be returned
	experiments <- data.frame(treatments = treatments, timePoints = timePoints, data = data)
	return(experiments)
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
			yAxisNoIntervention <- append(yAxisNoIntervention, dataList$data[[i]])
			
			# Get data for when there is intervention on the target
			} else if (dataList$treatments[[i]] == intervention) {
				
				yAxisIntervention <- append(yAxisIntervention, dataList$data[[i]])
				}
		}

	yAxisMax <- range(0, yAxisNoIntervention, yAxisIntervention)
	plot(xAxis, yAxisNoIntervention, xlab="Time Points", ylab=paste(gsub("DV.", "", childAntibody), "Data", sep=" "),
		type = "b", col = "blue", ylim = yAxisMax)
	lines(xAxis, yAxisIntervention, type = "b", pch = 22, col="red", lty = 2)
	legend("bottomright", c("No Intervention", "Intervention"),
		cex=0.8, col = c("blue", "red"), pch=21:22, lty=1:2, lwd=2, bty="n");
	title(main=paste(gsub("DV.", "", childAntibody), "Expression Data", sep=" "), col.main="red", font.main=4)
	} # End of dataPlotter function

dataList <- midasReader(midasFilePath, childAntibody)
dataPlotter(dataList, noIntervention, intervention)
