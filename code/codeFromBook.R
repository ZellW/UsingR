headtail(michelson)#headtail is specific to UsingR

wts  <- kid.weights[,2]

fivenum(wts)
#returns five number summary (minimum, lower-hinge, median, upper-hinge, maximum) for the input data.
range(wts)
diff(range(wts)) #
var(wts)
sd(wts)
z_scores <- scale(wts)[,1] #provides the z scores

#Example:  Professor grades by z scores.  A student with a z score >1.28 gets an A.  What grades is needed to get an A givent the score vector below?

grades <- c(54, 50, 79, 79, 51, 69, 55, 62, 100, 80)
answer <- mean(grades) + 1.28 * sd(grades) #ans is 88.9
iqr <- IQR(grades) #Inter Quartile Range - the distance or range between the Quartile 1 and Quartile 3
qqnorm(wts)
stripchart(wts) #next looks better
dotplot(wts) #gglpot can produce nice ones too http://www.r-bloggers.com/summarising-data-using-dot-plots/
boxplot(wts, kid.weights$height, horizontal=TRUE)

boxplot(Speed ~ Expt, data=michelson)#Model formula example. Speed is broken into groups defined my Expt (Speed in theis case is the y value) p 98 UsingR Book
#It shows how Speed (y) varies across groups (x) defined by Experiment
boxplot(Speed ~ Expt, data=michelson, subset=Expt %in% 3:4)
#subset displays only experiments 3 and 4 in the dataset - filtering example
#y ~ x:  y is the response (dependent) variable; x is the predictor (independent) variable
plot(fat$wrist, fat$neck) #is the ame as
plot(neck ~ wrist, fat)

plot(Speed ~ Expt, data=michelson)#compare to below
plot(summary(Speed ~ Expt, data=michelson))#Look at summary(Speed ~ Expt, data=michelson) and summary(height ~ weight, data=kid.weights) - interesting

#Paired Data (Numerical)P 102
cor(fat$wrist, fat$neck) #strongly correlated
cor(fat$wrist, fat$height) #mildly correlated
cor(fat$age, fat$ankle) #uncorrelated 

cor(ToothGrowth$dose, ToothGrowth$len) #strong correlation:  0.8026913
l <- split(ToothGrowth$len, ToothGrowth$dose)
groupmeans <- c(mean (l[[1]]), mean(l[[2]]), mean(l[[3]]))
cor(c(0.5, 1, 2), groupmeans)
#correlations formed from averages are typically closer to 1 or -1 than when all the data is considered individually  see page 11 of UsingR

cor(SAT$salary, SAT$total)#no correlation - is is Negative
#perc is the percentage of stidents that took the SAT - let's evaluate the effect of this:
plot(total ~ salary, data=SAT)# no correlation evident
points(total ~ salary, SAT, subset = perc<10, pch=15)#squares
points(total ~ salary, SAT, subset = perc>40, pch=16)#solid
#correlation for each subgroup appears positive
cor10 <- SAT %>% filter(perc<10) %>% select(salary, total)
corMiddle <- SAT %>% filter(perc>10, perc<40) %>% select(salary, total)
cor40 <- SAT %>% filter(perc>40) %>% select(salary, total)
c(less=cor(cor10$total, cor10$salary), middle=cor(corMiddle$total, corMiddle$salary), cor(cor40$total, cor40$salary))#returns all positive
#Correlations for all subgroups are positive when overall is negative.  
#Simpson's Paradox (when a trend for subgroups changes when data is aggregated.

#Linear Regression
residual <- lm(maxrate ~ age, data=heartrate)#to determine max heart rate by age, tradition is to 220- age (a line slope -1 intercept 220)
# residuals = observed - expected - sum of all residuals theortically add to 0
plot(maxrate ~ age, heartrate)
abline(residual)
#can use the predict function to calculate results from a lm:
predict(residual, data.frame(age=c(30, 40)))
plot(residual)#pretty cool -s ee all teh plots that are created

#Paired Bivariate Data:  categorical data (P 132)
#Usually involves counts of 2 categories - The distribution of an individual variable is called marginal distribution
#Sometimes you need to build the data - matrix is useful
seatbelts <- matrix(c(56, 2, 8, 16), nrow=2)
dimnames(seatbelts) <- list(parent=c("buckled", "unbuckled"), child=c("buckled", "unbuckled"))
#2 way tables can also be made with the table function
gradeTable <- table(grades$prev, grades$grade)
margin.table(seatbelts, margin=1)# margin 1 is rows; 2 for columns - this adds the row numbers "
margin.table(gradeTable, margin=2)
addmargins(seatbelts)

barplot(seatbelts, xlab="Parent", main="Child Seat Belt Usage") #stacked bars
barplot(seatbelts, xlab="Parent", main="Child Seat Belt Usage", beside=TRUE)
#to switch the data and use the rows to be the primary distribution, transpose:
t(barplot(seatbelts, xlab="Child"))
#to use proportions rather than counts for the y axis, use the prop.table
barplot(prop.table(seatbelts), xlab="Parent", main="Child Seat Belt Usage") #stacked bars

#mosaic plots are useful for 2 or more categorical variables (p 141)
#Using the Titanic data, it is a table that needs to be converted to a data frame (therwise all cols are chr type)
titanic <-  as.data.frame(Titanic)#includes 4 factors and 1 column of numbers for Suriviors
mosaicplot(xtabs(Freq ~ Sex, titanic))
mosaicplot(xtabs(Freq ~ Sex + Survived, titanic))
mosaicplot(xtabs(Freq ~ Sex + Survived + Class, titanic))
mosaicplot(xtabs(Freq ~ Class +Sex + Survived, titanic))#Order makes a difference!
mosaicplot(xtabs(Freq ~  Survived + Class, titanic))
#You can make ordered data out of factors
Survived <- rep(titanic$Survived, titanic$Freq)#makes a factor class of Yes and No for each of the 2201 passengers
Survived <- ordered(Survived) #ordered No then Yes
Class <- rep(titanic$Class, titanic$Freq)
Class <- ordered(Class)
#Ordered factors can be coerced into numeric data to allow fetting correlation:
cor(as.numeric(Survived), as.numeric(Class), method="kendall")
#negative answer is because of the ordering of the class with 1st being 1 and crew 4#kendall is used with cor with data 
#that can be ranked like ordered factors

#Multivariate Graphics p 204 UsingR
p4 <- ggplot(Cars93, aes(x=MPG.highway)) + stat_bin(binwidth=5)#binwidth has a defalt of length or range / 30
#Rather than display counts, you can display proportions - a density plot:
p5 <- ggplot(Cars93, aes(x=MPG.highway, y=..density..)) + stat_bin(binwidth=5)
#Explicit use of stat_bin is not required.  geom_histogram combines stat_bin and geom_bar
p6 <- ggplot(Cars93, aes(x=MPG.highway)) + geom_histogram(bin=5)
#Similarly, density plots can be made like this:
p7 <- ggplot(Cars93, aes(x=MPG.highway, y=..density..)) + geom_histogram(alpha=.5, binwidth=5) + geom_density()

p8 <- ggplot(Cars93, aes(x=Weight, y=MPG.highway)) + geom_point()+ geom_smooth()
p8 <- ggplot(Cars93, aes(x=Weight, y=MPG.highway)) + geom_point()+ geom_smooth(method="lm", se=FALSE)#removes the standard error shading from the plot above
#There are serveral methods.  Below we introduce a polynomial (p 207 UsingR):
p8 <- ggplot(Cars93, aes(x=Weight, y=MPG.highway)) + geom_point()+ geom_smooth(method="lm",formula=y ~ poly(x,2), se=FALSE)
#ggplot and faceting - see p 207 UsingR
#ggplot faceting uses the model formula syntax. .~f or f~. facets grouped by rows or columns, respectfully
p9 <- ggplot(PearsonLee, aes(y=child, x=parent)) + geom_point(alpha=.5) + 
     geom_smooth(method="loess") + facet_grid(par ~ chl)

Cars93_1 <- mutate(Cars93, price=cut(Price,c(0, 15, 30, 75),labels=c("Cheap", "Moderate", "Expensive")))
p10 <- ggplot(Cars93_1, aes(x=Weight, y=MPG.highway)) + geom_point(cex=3) + geom_smooth(method="lm", se=FALSE) +
     facet_grid(~ price)
#Margins - not sure what this really means (p 209 UsingR)
p11<- ggplot(PearsonLee, aes(y=child, x=parent)) + geom_point(alpha=.5) + 
     geom_smooth(method="loess") + facet_grid(par ~ chl, margins="chl")
#Facet Wrap is used when faceting by a factor:
p12 <- ggplot(morley, aes(x=Speed)) + geom_histogram(binwith=50) + facet_wrap(~Expt)#lists each of the 5 experiements
