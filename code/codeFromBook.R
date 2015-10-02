#See http://www.r-bloggers.com/?s=Using+R+for+Introductory+Statistics
#http://www.inside-r.org/packages/cran/UsingR/docs/getAnswer
#http://rpackages.ianhowson.com/cran/UsingR/
#
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

#Populations - Sample (P 219 UsingR)
sample1 <- sample(0:1, size=10, replace=TRUE)# toss of a coin.  default is to sample without replacement
sample2 <- sample(1:6, size=10, replace=TRUE)#roll of die 10 times
sample3 <- sample(1:6, size=10, replace=TRUE) + sample(1:6, size=10, replace=TRUE)
sample4 <-  sample(rep(0:1), times=c(3200, 6800), size=10, replace=TRUE) #rep(0:1, times=c(3200, 6800)) gives 3200 0's and 6800 1's in order

#Distribution Families (P 222 UsingR)
dunif(x=1, min=0, max=3) # d provides answers like what is the probability of x=5 - not 5 or less like below
punif(q=2, min=0, max=3) # p cumulative probability function for binomial distribution pbinom (gives answers like 5 or less)
qunif(p=1/2, min=0, max=3) # q returns the quartile - obviously 1/2 of 3 is 1.5
runif(n=1, min=0, max=3) # r returns a random value in [0,3]
runif(10, min=0, max=1:5) #returns 10 random values between 0 and 5

#See http://www.cyclismo.org/tutorial/R/probability.html
# dnorm. Given a set of values it returns the height of the probability distribution at each point. 
# If you only give the points it assumes you want to use a mean of zero and standard deviation of one. 
# 
# pnorm. Given a number or a list it computes the probability that a normally distributed random 
# number will be less than that number. This function also goes by the rather ominous title of the 
# “Cumulative Distribution Function.” It accepts the same options as dnorm
#  you wish to find the probability that a number is larger than the given number you can use the lower.tail option
# 
# qnorm which is the inverse of pnorm. The idea behind qnorm is that you give it a probability, and 
# it returns the number whose cumulative distribution matches the probability. For example, if you 
# have a normally distributed random variable with mean zero and standard deviation one, then 
# if you give the function a probability it returns the associated Z-score

x <- runif(100) #100 random values from uniform distribution
d <- density(x) #note the outpunt if d is a class of density! This is what is returned:

# Call:
#      density.default(x = x)
# 
# Data: x (100 obs.);	Bandwidth 'bw' = 0.09965
# 
# x                 y           
# Min.   :-0.2925   Min.   :0.001617  
# 1st Qu.: 0.1047   1st Qu.:0.160913  
# Median : 0.5019   Median :0.706437  
# Mean   : 0.5019   Mean   :0.628747  
# 3rd Qu.: 0.8991   3rd Qu.:1.033379  
# Max.   : 1.2962   Max.   :1.257553

curve(dunif, -.1, 1.1, ylim=c(0,max(d$y,1))) #plots the uniform distribution
curve(dunif, min(d$x), max(d$x), ylim=c(0,max(d$y,1))) #my version ;)
lines(d, lty=2) # adds the density estimate - the curve would follow the uniform distribution with an increase of random samples
rug(x) #adds a rug - 1-d plot - to the plot

#Bernoulli Random Values - only 2 possible values - 0 or 1
# R generates samples from this distribution using sample:
n <- 10; p <- 1/4
sample(0:1, size=n, replace=TRUE, prob=c(1-p, p)) #Bernoulli distribution when the probability is 25%

#Binomial Random Values - counts the number of succsses in a Bernoulli trial
# See http://www.r-tutor.com/elementary-statistics/probability-distributions/binomial-distribution AND the RMD in this project
#Toss a coin 10 times
dbinom(5, size=10, prob=1/2) #result is there is a 24.6% probability that 5 heads will be counted
#probability that there are 6 or fewer:
sum(dbinom(0:6, size=10, prob=1/2))
#or
pbinom(6, size=10, p=1/2)

#Normal Distribution (P227 UsingR)
# 2 normal distributions have the same probabilities (area under hte curve to the left) at a given sd:
pnorm(1, mean=0, sd=1) # 1 is the value where the curve is centered on the x axis - the answer is the z-score
pnorm(4.5, mean=4, sd=.5) # 4.5 is the value where the curve is centered on the x axis- the answer is the same z-score (84.1% o fthe area is to the left)

qnorm(c(.25, .5, .75))#Gets the z-scores for the 2nd and 3rd quartiles - the IQR.  IQR ~ 1.35sigma

pnorm(1)-pnorm(-1) #tells ue that ~68% of the area is between 1 sd of the mean - for any normal distribution
pnorm(2)-pnorm(-2) #tells ue that ~95% of the area is between 2 sd of the mean - for any normal distribution
pnorm(3)-pnorm(-3) #tells ue that ~99.7% of the area is between 3 sd of the mean - for any normal distribution

#What percent of men are at least 6 foot given mean=70.2 and sd=2.89?
mu <- 70.2; sigma <- 2.89
pnorm(6*12, mean=mu, sd=sigma) #73% men are 6 foot or shorter

#How tall is the tallest man?
#There are about 3.5 billion males.  So 1 in 3.5B:
prob <- 1-1/3500000000
qnorm(prob, mu, sigma)/12 #about 7.34 feet tall - not very accurate with only 2 paremeters but still a good exercise

#Uniform Distribution (p 231 UsingR)
res <-  runif(50, min=0, max=10)#Get random vaiables from uniform distribution
par(fig=c(0,1,0,.35))#Sets up the the plot placement using the lower 35% of the diagram
boxplot(res, horizontal=TRUE, bty="n", xlab="Uniform Sample")#the type of box to be drawn around the legend. The allowed values are "o" (the default) and "n".
par(fig=c(0,1,.25,1), new=TRUE)# fig= setting uses top- 75% of figure
hist(res, prob=TRUE, main="", col=gray(.9))
lines(density(res), lty=2)
curve(dunif(x, min=0, max=10), lwd=2, add=TRUE)
rug(res)

#Exponential Distribution
#See http://www.r-bloggers.com/using-r-for-introductory-statistics-chapter-5-probability-distributions/
samples <- rexp(100, rate=1/5)
par(fig=c(0,1,0,0.35))
boxplot(samples, horizontal=T, bty="n", xlab="log-normal distribution")
par(fig=c(0,1,0.25,1), new=T)
s <- seq(0,max(samples),0.1)
d <- dexp(s, rate=1/5)
hist(samples, prob=T, main="", col=gray(0.9), ylim=c(0,max(d)))
lines(density(samples), lty=2)
curve(dexp(x, rate=1/5), lwd=2, add=T)
rug(samples)


#Lognormal Distribution (p 233 UsingR)
#A heaveliy skewed continuous distribution on positive number.  Can describe income or survival distributions as an example
#Use lnorm with meanlog and sdlog as pararmeters

samples <- rlnorm(100, meanlog=0, sdlog=1)
par(fig=c(0,1,0,0.35))
boxplot(samples, horizontal=T, bty="n", xlab="log-normal distribution")
par(fig=c(0,1,0.25,1), new=T)
s <- seq(0,max(samples),0.1)#numerical vector from 0 to the max of samples - 10.6 incremented by .1
hist(samples, prob=T, main="", col=gray(0.9), ylim=c(0,max(d)))
lines(density(samples), lty=2)
curve(dlnorm(x, meanlog=0, sdlog=1), lwd=2, add=T)
rug(samples)

#t, f and chisq distributions containing 95% of the area with degrees of freedom = 10

qt(c(.025, .975), df=10)
qf(c(.025, .975), df1=10, df2=5)#f distribution requires 2 df values
qchisq(c(.025, .975), df=1)

#Central Limit Theorem (p 236 UsingR)
#Grocery checkoutline.  Mean service mu = 1; std dev is 1 minute. What is prob that ave checkout times be .9 minutes or less?
pnorm(.9, mean=1, sd=1/sqrt(20))#32.7%

#Not in book, found online at http://rpubs.com/aousabdo/BOC1
x <- seq(0,50,1)
y <- dbinom(x,50,0.2)#recalll vector, number of trials, probability
plot(x,y, col="blue", "l", ylab="Binomial Density")

#Statistical Inference - Simulations (p 244 UsingR)
#Use replicate
#M is the number of times you run the sample of 16 - testing 10 times with sample size 16
func1 <- function(M, n, mu, sd) {replicate(M, mean(rnorm(n, mean=mu, sd=sigma)))}
func1(10, 16, 100, 16)# produces M (10) realizations of sample (size=16) mean for normal population
#Find z-score of a sample
zstat <- function(x, mu, sigma){
     (mean(x)-mu)/(sigma/sqrt(length(x)))
}
#Replicate the z-score - replicate function take a value n and an expression
M <- 2000; n <- 7; mu <- 100; sigma <- 16
res <- replicate(M, {
     x <- rnorm(n, mu, sigma)
     zstat(x, mu, sigma)
})
hist(res)
#t distribution
tstat <- function(x, mu){
     (mean(x)-mu) / (sd(x) /sqrt(length(x)))}
mu <- 0; sigma <- 1; M <- 750; n <- 4
res2 <- replicate(M, tstat(rnorm(n, mu, sigma), mu))
boxplot(res2)
hist(res2)

#Mean vs Media - good plottintg exercise
M <- 1000; n <- 35
res_mean <- replicate(M, mean(rnorm(n)))
res_median <- replicate(M, median(rnorm(n)))
boxplot(list("Sample Mean"=res_mean, "Sample Median"=res_median), main="Normal Population")
#Note the median have greater variability

#Estimating Probabilities (p 250 UsingR)
#Example - Lottery - N diff numbers, player selects k number, j numbers are selected at random; player wins if 1 or more matched
N <- 80; k <- 20; j <- 10#out of 80 numbers, player selects 20numbers and hopes it matches one of the 10 numbers the computer selects from the 80
x <- sample(1:N, j, replace=FALSE)
#Let's see if you won - were number 1 through 20 one of the 10 selected from 80?
sum(x %in% 1:k)
#simulate  1000 times
M <-  10000
res3 <- replicate(M, {
     x <- sample(1:N, j, replace=FALSE)
     sum(x %in% 1:k)
})
#How many times were 0 or 1 matches found in a population on 10,0000 simulations
c("Zero Matches"=sum(res3==0)/length(res3), "One Match"=sum(res3==1)/length(res3))
#Liklihood of 5 or more matches?
sum(res3>=5)/length(res3)

#Significance Test (p 252 UsingR)
#Does honey improve performance? 7 participants; 3 in control groiup and 4 in treatement.  The data is:
testCtrl <- c(23, 33, 40); testTreatment <- c(19, 22, 25, 26)
the_data <- stack(list(ctrl=testCtrl, treatment=testTreatment))#Stacking vectors concatenates multiple vectors into 
#a single vector along with a factor indicating where each observation originated. Unstacking reverses this operation
aggregate(values~ind, the_data, mean)# or
the_data %>% group_by(ind) %>% summarise(TheMean=mean(values))
#Perhaps the results (the mean with honey was 9 higher) were do to who was chosen.  Let's simulate
cmbs <- combn(7,3) #Generate all combinations of the elements of x taken m at a time. If x is a positive integer, 
#returns all combinations of the elements of seq(x) taken m at a time. 
#So it creates all the combinations of the 7 participants 3 at a time = 35 possible unique combinations - order does not matter
#To know how many combinations rather than developing all the data is:
cmbs_count <- choose(7,3)
#The first case in the matrix would be what we already observed - lets demo this:
ind <- cmbs[,2]#is the first set of 3 so the referenced values are still 23, 33, 40 
obs <- mean(the_data$value[ind] - mean(the_data$value[-ind]))#selects the first 3 and then selects everything else left
#Let's simulate:
res4 <- apply(cmbs, 2, function(ind){
     mean(the_data$values[ind]) - mean(the_data$values[-ind])
     })#the values represent the radomization distribution for the difference in group means - how extreme is the value 9?
sum(res4>=obs)#3
sum(res4>=obs)/length(res4)#8.6%
#3 of 35 wil be equal or larger than 9. Seems unlikly but not impossible that honey improves performance.
hist(res4)#just for fun

#p 254 UsingR
#Does caffiene make you jittery?  20 classmates, 10 in each team.  Every other is service caffinated coffee.  Count finger taps.  The results:
caf <- c(245, 246, 246, 248, 248, 248, 250, 250, 250, 252)
no_caf <- c(242, 242, 242, 244, 244, 245, 246, 247, 248, 248)
the_data2 <- stack(list(Caffiene=caf, No_Caffiene=no_caf))
obs2 <- mean(caf)- mean(no_caf)#3.5
#How unusual is thhis if all possible random assingment of the 20 into 2 groups?
tmp2 <- choose(20,10)#184,756 - a bit much.  Consider if the set was larger:
tmp3 <- choose(60,30)#HUGE number - too much even for computers.  Alternative:  simulate from the radominzation distribution
sample(1:20, 10)#Randomly select 10
#Now simulate this
res5 <- replicate(2000, {
     ind <- sample(1:20, 10, replace=FALSE)
     mean(the_data2$values[ind]) - mean(the_data2$values[-ind])
})
sum(res5>obs2)/length(res5)#.0025 ~ 0.25% - therefore strong indication caffiene increases finger tapping

#How well does a sample statistic estimate a parameter? (Above we wanted to know if a treatment induced an effect)
#Bootstrap
#Sample data  on diff between charge and Medicare payment fopr a diagnoses:
data (Medicare)
diabetes <- subset(Medicare, subset=DRG.Definition == "638 - DIABETES W CC") #or
diabetes <- filter(Medicare, DRG.Definition == "638 - DIABETES W CC")
priceGap <- with(diabetes, Average.Covered.Charges - Average.Total.Payments)
range(priceGap)#Very large range
#Goal is to find asn interval where we have conficence the unknown population mean resides
xbar <- mean(priceGap)
#A single random sample:
xstar <- sample(priceGap, length(priceGap),replace=TRUE)
#Bootstrap uses these values as if they were actually from the population.  Example:
mean(xstar-xbar)
#So lets replciate.
M <- 2000; 
res6 <- replicate(M,{
     xstar <- sample(priceGap, length(priceGap), replace=TRUE)
     mean(xstar)- xbar
})
#Let's npw use this to get 95% conficence:
alpha <- .05
xbar + quantile(res6, c(alpha/2, 1-alpha/2))

#Bayesian Analysis -p 259 UsingR - SKIPPED for COURSERA CLASS
#
#Confidence Intervals - p 262 UsingR
#Proportions######################################
#Random sample 80 students. 80 had moderate exercise 5 or more times a week. What is 90% confidence level for the proportion of all 1,812 students?
#Manual Method:
x <- 80; n <- 125
proportion <- x/n
alpha <- 1-0.90
zstar <- pnorm(1-alpha/2)#provides the quantile density function - why quantile - not sure - pnorm gives the same results.  
#See plot p 267 to understand 1-alpha
SE <- sqrt(proportion*(1-proportion)/n)
MOE <- zstar*SE #Margin Of Error
proportion + c(1,-1)* MOE#0.7106177 0.5693823
#Compare using R formula
prop.test(x,n, conf.level = .9)#same result
prop.test(x,n, conf.level = .9)$conf.int#to get just the part we are interest in

#Survey says 57%, n=1000, MOE=3%
#Not that hard.  We know MOE=z*SE; z=MOE/SE; z=.03/SE
z <- .03/sqrt(.57*(1-.57)/1000)
alpha <- 2*pnorm(-z)# 
(1-alpha)*100#implied 95%

#Confidence Internvals - Population Mean p 271 UsingR###########################
#t test
#NOTE:  SE(proportions) = sqrt((P*(1-P)/n); SE(mean) = sd/sqrt(n)
#Class of 30 average height 66, sd=4. What is 80% conf level?
xbar <- 66; sd <- 4; n <- 4; alpha <- 1-.8
tstat <- qt(1-alpha/2, df=n-1)
SE <- sd/sqrt(n)
MOE <- tstat*SE
xbar+c(1,-1)* MOE

