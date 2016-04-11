# Add a slash to the end of a directory listing, if required.
add.slash = function(path) {
  if(substring(path, nchar(path), nchar(path)) != "/") {
    path = paste(path, "/", sep = "")
  } # if(substring(path, nchar(path), nchar(path)) != "/"))
  return(path)
} # add.slash()
