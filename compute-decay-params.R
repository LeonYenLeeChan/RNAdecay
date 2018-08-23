#!/usr/bin/Rscript --vanilla --slave

## computes half-life and efficiency parameters for a decay dataset

## usage: Rscript compute_decay_params.R  input

## takes an input file with the following tab delimited fields:
##	filenames	the names of the fpkm_tracking files (cufflinks output) that contain the fpkm data for each timepoint - list the files in increasing time order
##	timepoints	the times at which the samples were collected at
##	spike		the name of the gene_id of the spike RNA in the gff file - standard ones are hSRP1 and xlRCC1
##	keep_LOWDATA	a logical that includes data from cufflinks LOWDATA entries
##	cutoff_file	a .fpkm_tacking file to use to filter out low abundance transcripts
##	FPKM_cutoff	the value below which the gene is considered low abundance in the cutoff_file
##	Rsq_cutoff	the value below which the fit is considered bad

## outputs 4 (or 5) files:
##	decaydata.raw			a tab delimited txt file of just the raw data
##	decaydata.spike.norm		a tab delimited txt file of the raw data divided by the spike values
##	decaydata.norm			a tab delimited txt file of the spike normalized data renormalized such that the first timepoint is = 1
##	decaydata.fit			a tab delimited txt file which contains the fitted halflife, the fitted efficiency and the Rsq


## load libraries
require("MASS")
require("minpack.lm")
require("rgl")
require("robustbase")
require("qpcR")

## get data from input file and assign to variables
cmdArgs <- commandArgs(trailingOnly = TRUE) ##gets the command line arguements
input_nlines <- as.numeric(strsplit(system(paste("wc -l ", cmdArgs[1], sep = ""), intern = TRUE), " ")[[1]][1])

for (i in 1:input_nlines)
	{
		ith_line <- scan(cmdArgs[1], what = "", skip = i-1, nlines = 1)
		if (ith_line[1] == "filenames") filenames <- ith_line[2:length(ith_line)]
		if (ith_line[1] == "timepoints") timepoints <- as.numeric(ith_line[2:length(ith_line)])
		if (ith_line[1] == "spike") spike <- ith_line[2]
		if (ith_line[1] == "keep_LOWDATA") keep_LOWDATA <- as.logical(ith_line[2])
		if (ith_line[1] == "cutoff_file") cutoff_file <- ith_line[2]
		if (ith_line[1] == "FPKM_cutoff") FPKM_cutoff <- as.numeric(ith_line[2]) else FPKM_cutoff <- 2
		if (ith_line[1] == "Rsq_cutoff") Rsq_cutoff <- as.numeric(ith_line[2]) else Rsq_cutoff <- 0.8
	} #end for i

## reads each fpkm_tracking into a dataframe and stores it as an element of the "rawdata" list
rawdata <- list()
for (i in 1:length(filenames))		rawdata[[i]] <- read.table(filenames[i], header=TRUE, sep="", fill=TRUE)

## remove entries from "rawdata" that do not have a gene_id
for (i in 1:length(rawdata))		rawdata[[i]] <- rawdata[[i]][rawdata[[i]]$gene_id != "-",]

rawdata_noQC <- rawdata ## generate list of the raw data where LOWDATA and HIDATA are not taken into account
  
## replace fpkm values with placeholders for FPKM_status with HIDATA or LOWDATA
for (i in 1:length(rawdata)) 
	{
		if (!keep_LOWDATA) rawdata[[i]]$FPKM[rawdata[[i]]$FPKM_status=="LOWDATA"] <- -1
		rawdata[[i]]$FPKM[rawdata[[i]]$FPKM_status=="HIDATA"] <- -2
	}

## read in cutoff file and clean up data
cutoff_data <- read.table(cutoff_file, header = TRUE, sep = "", fill = TRUE)
cutoff_data <- cutoff_data[cutoff_data$gene_id != "-", ]
cutoff_data$FPKM[cutoff_data$FPKM == "LOWDATA"] <- -1
cutoff_data$FPKM[cutoff_data$FPKM == "HIDATA"] <- -2

## merge fpkm data into one data frame
## returns "decaydata" and "decaydata_noQC" dataframes which has 2 + the number of files columns, the first two being gene_id and gene_name and the last columns being the fpkm values
FPKMnames <- paste("T",1:length(rawdata),sep="")

for (i in 1:length(rawdata))
  {
    colnames(rawdata[[i]])[10] <- FPKMnames[[i]]
    colnames(rawdata_noQC[[i]])[10] <- FPKMnames[[i]]
    if (i == 1) {
	decaydata <- rawdata[[i]][, c(4, 5, 10)]
	decaydata_noQC <- rawdata_noQC[[i]][, c(4, 5, 10)]
		}
    else
      {
	decaydata <- merge(decaydata, rawdata[[i]][, c(4,10)])
	decaydata_noQC <- merge(decaydata_noQC, rawdata_noQC[[i]][, c(4,10)])
      }
  }

## remove temporary data frames
rm (filenames, rawdata)

## writes decaydata.raw
write.table(decaydata,"decaydata.raw",sep="\t",row.names=FALSE, quote = FALSE)
write.table(decaydata_noQC, "decaydata_noQC.raw", sep = "\t", row.names = FALSE, quote = FALSE)

#make an index of bad values (fpkm equals -1 or -2) such that these values do not get computed on
    bad_matrix <- decaydata[,3:(2+length(timepoints))] == -1 | decaydata[,3:(2+length(timepoints))] == -2
    bad_index <- logical(0)
    for (i in 1:length(bad_matrix[,1])) bad_index <- c(bad_index, TRUE %in% as.vector(bad_matrix[i, ]))                                                                                       

## normalize data according to norm_method
    decaydata_spike_norm <- decaydata
    decaydata_norm <- decaydata
    spike_vector <- as.vector(decaydata[decaydata$gene_id == spike, 3:(2+length(timepoints))], mode = "numeric")
    spike_vector <- spike_vector/max(spike_vector)
    for (i in 1:dim(decaydata)[1])
      {
        if (!bad_index[i])
		{
			##divide fpkm values by spike values
			decaydata_spike_norm[i, 3:(2+length(FPKMnames))] <- as.vector(decaydata[i, 3:(2+length(timepoints))], mode = "numeric")/spike_vector
       			##divide each fpkm value by T1 value
			decaydata_norm[i, 3:(2+length(FPKMnames))] <- as.vector(decaydata_spike_norm[i, 3:(2+length(timepoints))], mode = "numeric")/decaydata_spike_norm[i,3]
		} #end if !bad_index[i]
      } ##end for (i in 1:dim(decaydata)[1])

## write spike snormalized data to file
write.table(decaydata_spike_norm,"decaydata.spike.norm",sep="\t",row.names=FALSE, quote = FALSE)

## write normalized data to file
write.table(decaydata_norm,"decaydata.norm",sep="\t",row.names=FALSE, quote = FALSE)


##start here!!
## compute fitted values for decay data using exponential decay model and the  plateau-exponential decay model
decaydata_fitted <- decaydata_norm[,1:2]
decaydata_fitted[,3:20] <- NA
colnames(decaydata_fitted)[3:20] <- c("Th_final", "eff_final", "Rsq_final", "Th_fit", "eff_fit", "Rsq", "lim_Th_fit", "lim_eff_fit", "lim_Rsq", "data_has_0s", "FPKM_LT_cutoff",  "data_has_nonfinites", "data_has_bad_data", "decay_below_exp_threshold", "second_timepoint_GT_95xfirst_timepoint", "exp_plateau_fit_successful", "limited_exp_plateau_fit_successful", "Rsq_LT_cutoff") 
 
for (i in 1:dim(decaydata_norm)[1])
    {
    	## make notes but compute values all the same values
    	normdata <- as.vector(decaydata_norm[i, 3:(2+length(timepoints))], mode = "numeric")
    	ith_timepoints <- timepoints
    	ith_try_lim <- FALSE
    	if (i%%500 == 0) print(paste(i, " genes analyzed...", sep = ""))
    	if (i == dim(decaydata_norm)[1]) print(paste(i, " total genes analyzed", sep = ""))
	##test FPKM value in cutoff file is less than cutoff_FPKM value
        if (cutoff_data$FPKM[cutoff_data$gene_id == decaydata_fitted$gene_id[i]] < FPKM_cutoff) decaydata_fitted$FPKM_LT_cutoff[i] <- TRUE else decaydata_fitted$FPKM_LT_cutoff[i] <- FALSE
	##test if FPKM vector contains zero values
	if (0 %in% normdata) decaydata_fitted$data_has_0s[i] <- TRUE else decaydata_fitted$data_has_0s[i] <- FALSE 
        ##test if FPKM vector contains non-finites
	if (!all(is.finite(normdata))) decaydata_fitted$data_has_nonfinites[i] <- TRUE else decaydata_fitted$data_has_nonfinites[i] <- FALSE
	##test if FPKM vector contains bad data
	if (bad_index[i]) decaydata_fitted$data_has_bad_data[i] <- TRUE else decaydata_fitted$data_has_bad_data[i] <- FALSE
	##try fits if no nonfinites and no bad data
	if (!decaydata_fitted$FPKM_LT_cutoff[i] && !decaydata_fitted$data_has_nonfinites[i] && !decaydata_fitted$data_has_bad_data[i])
		{
			if (normdata[2] > 0.95*normdata[1]) ## note and trim first timepoint if FPKM at second timepoint is > 0.95 of FPKM at first timepoint
        		        {
                        		##renormalize to second timepoint
                        		ith_timepoints <- (timepoints - timepoints[2])[2:length(timepoints)]
                        		normdata <- (normdata/normdata[2])[2:length(normdata)]
                        		decaydata_fitted$second_timepoint_GT_95xfirst_timepoint[i] <- TRUE
                		} else decaydata_fitted$second_timepoint_GT_95xfirst_timepoint[i] <- FALSE
			##test if FPKM vector has the form of (1,0,0, ...)
          		if(sum(normdata[2:length(normdata)]) == 0) decaydata_fitted$decay_below_exp_threshold[i] <- TRUE else decaydata_fitted$decay_below_exp_threshold[i] <- FALSE
                   	if(!decaydata_fitted$decay_below_exp_threshold[i])
				{
					decay_dataframe <- data.frame(normdata, ith_timepoints)
					## plateau-exponential model
      					##fit without limits
      					fit <- try(nls(normdata~eff*2^(-ith_timepoints/Th) + (1-eff), data = decay_dataframe, start = list(eff = 1, Th = 1)), silent = TRUE)
      					##fit with limits
      					limited_fit <- try(nls(normdata~eff*2^(-ith_timepoints/Th) + (1-eff), data = decay_dataframe, start = list(eff = 1, Th = 1), algorithm = "port", lower = list(eff = 0,Th = 0), upper = list(eff = 1, Th = 1000000)), silent = TRUE)
      					decaydata_fitted$exp_plateau_fit_successful[i] <- FALSE
					decaydata_fitted$limited_exp_plateau_fit_successful[i] <- FALSE
					##write fit data
      					if (class(fit) != "try-error")
        					{
							fit_summary <- try(summary(fit), silent = TRUE)
							if (class(fit_summary) != "try-error")
								{
         								decaydata_fitted$Th_fit[i] <- fit_summary$coefficients[2,1]
         								decaydata_fitted$eff_fit[i] <- fit_summary$coefficients[1,1]
         								decaydata_fitted$Rsq[i] <- Rsq(fit)
 									decaydata_fitted$exp_plateau_fit_successful[i] <- TRUE 
								}
	       					}
     					##write limited fit data
      					if (class(limited_fit) != "try-error")
        					{
							limited_fit_summary <- try(summary(limited_fit), silent = TRUE)
							if (class(limited_fit_summary) != "try-error")
								{
         								decaydata_fitted$lim_Th_fit[i] <- limited_fit_summary$coefficients[2,1]
         								decaydata_fitted$lim_eff_fit[i] <- limited_fit_summary$coefficients[1,1]
         								decaydata_fitted$lim_Rsq[i] <- Rsq(limited_fit)
       									decaydata_fitted$limited_exp_plateau_fit_successful[i] <- TRUE
								}
						}
					##determine final Th, eff and Rsq values
					if (decaydata_fitted$exp_plateau_fit_successful[i])
						{
							if (decaydata_fitted$Rsq[i] > Rsq_cutoff)
								{
									if (decaydata_fitted$eff_fit[i] < 1)
										{
											decaydata_fitted$Th_final[i] <- decaydata_fitted$Th_fit[i]
											decaydata_fitted$eff_final[i] <- decaydata_fitted$eff_fit[i]
											decaydata_fitted$Rsq_final[i] <- decaydata_fitted$Rsq[i]
											decaydata_fitted$Rsq_LT_cutoff[i] <- FALSE
										}
									else ith_try_lim <- TRUE 
								}
							else ith_try_lim <- TRUE
						}
					else ith_try_lim <- TRUE
					if (ith_try_lim)
						{
							if (decaydata_fitted$limited_exp_plateau_fit_successful[i])
								{
									if (decaydata_fitted$lim_Rsq[i] > Rsq_cutoff)
										{
											decaydata_fitted$Th_final[i] <- decaydata_fitted$lim_Th_fit[i]
											decaydata_fitted$eff_final[i] <- decaydata_fitted$lim_eff_fit[i]
											decaydata_fitted$Rsq_final[i] <- decaydata_fitted$lim_Rsq[i]
											decaydata_fitted$Rsq_LT_cutoff[i] <- FALSE
										}
									else decaydata_fitted$Rsq_LT_cutoff[i] <- TRUE
								}
						}
				}
		}
  } ## end for i loop

## write fitted values to file
write.table(decaydata_fitted, "decaydata.fit", sep = "\t", row.names = FALSE, quote = FALSE)
