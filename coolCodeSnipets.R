#nice way to identidy missing values:
apply(data,2,function(x) sum(is.na(x)))

beets <- c(41, 40, 41, 42, 44, 35, 41, 36, 47, 45)
no_beets <- c(51, 51, 50, 42, 40, 31, 43, 45)
#Stack data: p 99 of UsingR
b <- list("beets" = beets, "no beets"=no_beets)
b
stacked <- stack(b) #results in a data frame
plot(values ~ ind, data=stacked)

# Split function
speeds <- split(michelson$Speed, michelson$Expt)# date followed by factor to split by
speeds <- stack(speeds) #puts it back together
