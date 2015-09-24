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

#Correlation
cor(fat$wrist, fat$neck) #strongly correlated
cor(fat$wrist, fat$height) #mildly correlated
cor(fat$age, fat$ankle) #uncorrelated 

cor(ToothGrowth$dose, ToothGrowth$len) #strong correlation:  0.8026913
l <- split(ToothGrowth$len, ToothGrowth$dose)
groupmeans <- c(mean (l[[1]]), mean(l[[2]]), mean(l[[3]]))
cor(c(0.5, 1, 2), groupmeans)
#correlations formed from averages are typically closer to 1 or -1 than when all the data is considered individually  see page 11 of UsingR
