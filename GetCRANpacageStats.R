types <- c("source", "win.binary", 
           "mac.binary", "mac.binary.mavericks")

CRANmirror <- "http://cran.revolutionanalytics.com"

pdb <- lapply(types, function(x){
     cran <- contrib.url(repos = CRANmirror, 
                         type = x)
     available.packages(contriburl = cran, type = x)
})
names(pdb) <- types
str(pdb, max.level = 1)


# Number of available packages
sapply(pdb, nrow)

# Display first few packages for each type
pkgs <- sapply(pdb, rownames)
lapply(pkgs, head)

# Create list of differences between source and win.binary
setdiff(pkgs[["source"]], pkgs[["win.binary"]])
