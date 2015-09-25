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

seatbelts <- matrix(c(56, 2, 8, 16), nrow=2)
rownames(seatbelts)  <- c("buckled", "unbuckled")
colnames(seatbelts) <- c("buckled", "unbuckled")
#better way:
dimnames(seatbelts) <- list(parent=c("buckled", "unbuckled"), child=c("buckled", "unbuckled"))

gradeTable <- table(grades$prev, grades$grade)
margin.table(seatbelts, margin=1)# margin 1 is rows; 2 for columns - this adds the row numbers "
margin.table(gradeTable, margin=2)
addmargins(seatbelts)

#Automatically chage table data into a proportion by column or row:
table(grades$prev, grades$grade)
prop.table(table(grades$prev, grades$grade), margin=1)*100 #Devides each cell by the total count for the row where the cell is like 15/28 for cell (1,1)

xtabs(count ~ Whorls + Loops, Fingerprints)
plot(count ~ Whorls + Loops, Fingerprints) #did it just for fun
xtabs(~ Origin + Type, Cars93)#since nthing on left hand side or formula, it automatically tallies the data!
# the period - . - can be usedas shorthand notation for all variables in the data set not speficied on theleft hand side.
xtabs(count ~ ., Fingerprints)#since there are only 3 firleds in the dataset and count is specified then the 2 remaining fields (Whorls and Loops) 
#are assumed by the period shorthand



 