get.call.rate <-
function(path = ".", title = "", plot = F) {
  
  # Find a file with "DNAReport" in the filename.
  dnafile = dir(path = path, pattern = "DNAReport", full.names = T)

  # If not found, then quit.
  if(length(dnafile) == 0) {
    stop(paste("No file with 'DNAReport' in the filename was found in directory",
         path, ".  Please make sure that the DNAReport file is unzipped and",
         "in the specified directory."))
  } # if(length(dnafile) == 0)

  # If there is more than one call rate file, then quit.
  if(length(dnafile) > 1) {
    stop(paste("There is more than one file with DNAReport in the filename.",
         "Please place only one data set in each directory."))
  } # if(length(dnafile) > 1)

  # Read in the DNA file.  The current format needs us to skip two lines because
  # there is no comment character at the beginning of those lines.
  dnafile = read.csv(dnafile, skip = 2, stringsAsFactors = F)

  # Make sure that we have a "DNA_Name" column.
  if(length(grep("DNA_(Name|ID)", colnames(dnafile))) == 0) {
    stop(paste("There is no column called 'DNA_Name' in the DNAReport.",
         "Please label the column with the sample ID 'DNA_Name' in the",
         "DNAReport file and re-run the function."))
  } # if(length(grep("DNA_Name", colnames(dnafile))) == 0)

  # Make sure that we have a "Call_Rate" column.
  if(length(grep("Call_Rate", colnames(dnafile))) == 0) {
    stop(paste("There is no column called 'Call_Rate' in the DNAReport.",
         "Please label the column with the call rates 'Call_Rate' in the",
         "DNAReport file and re-run the function."))
  } # if(length(grep("DNA_Name", colnames(dnafile))) == 0)
  
  # Plot the call rate for all samples.
  if(plot) {
    png("call.rate.png", width = 800, height = 600)
    par(font = 2, font.axis = 2, las = 2)
    barplot(dnafile$Call_Rate, names.arg = dnafile$DNA_Name, ylim = c(0,1),
            main = paste("Call Rate:", title))
    abline(h = 0.9, col = 2)
    abline(h = 0.95, col = 2)
    dev.off()
  } # if(plot)

  # Return the call rates.
  retval = dnafile$Call_Rate
  names(retval) = dnafile$DNA_Name
  return(retval)

}
