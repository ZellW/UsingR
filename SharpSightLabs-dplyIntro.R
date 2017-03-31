#################################################
# DATA MANIPULATION WITH DPLYR
#
# How to ...
#   - filter rows:                  filter()
#   - select columns:               select()
#   - add new variables:            mutate()
#   - sort your data:               arrange()
#   - summarize your data:          summarize()
#   - group / aggregate your data:  group_by()
#
# (Release: 2016-08)
#
# sharpsightlabs.com
# © Copyright 2016 Sharp Sight Labs
# All rights reserved
#
#################################################


#----------------
# LOAD LIBRARIES
#----------------

library(dplyr)
library(ggplot2)


#-----------------
# DATA INSPECTION
#-----------------

# get first few records
head(diamonds)


# examine the levels of the "cut" variable
levels(diamonds$cut)


# create a frequency table of the "cut" variable
table(diamonds$cut)



############################################
# FILTER
# - filter() subsets your data and
#   keeps rows matching specific criteria
#
############################################


#---------------------------------------------------------
# "filter" rows of the diamonds dataset
#  and keep only records where the diamond cut is "Ideal"
#---------------------------------------------------------

# filter()
df.diamonds_ideal <- filter(diamonds, cut == 'Ideal')


# check first few rows
df.diamonds_ideal


# Check the new levels
levels(df.diamonds_ideal$cut)
table(df.diamonds_ideal$cut)


# Drop levels that no longer appear in the cut variable
df.diamonds_ideal <- droplevels(df.diamonds_ideal) 

head(df.diamonds_ideal)
levels(df.diamonds_ideal$cut)



#########################################
# SELECT
# - select() selects and keeps variables
#
#########################################


#-----------------------------------------------
# "select" specific variables from the diamonds
# dataset using select()
#-----------------------------------------------

# inspect first few rows
head(diamonds)


# select 1 columns: carat
df.diamonds_1variable <- select(diamonds, carat)

# inspect
head(df.diamonds_1variable)
names(df.diamonds_1variable)


# select 4 columns: carat, cut, color, price
df.diamonds_4variable <- select(diamonds, carat, cut, color, price)

# inspect
head(df.diamonds_4variable)




######################################################################
# MUTATE
# - mutate() creates a new variable in your dataframe
#    based on existing variables
# - i.e., creates a variable (which you give a name), and 
#   assigns a value; 
#   thus, the mutate() function takes name/value pairs as an argument
#
######################################################################

#-----------------------------------------------
# "mutate" the dataframe
# i.e., add a new variable to the diamonds data
#  using mutate()
#-----------------------------------------------

# inspect
head(diamonds)

# add a variable named "price per carat" to the data
df.diamonds_mutated <- mutate(diamonds, price_per_carat = price / carat)

# inspect
head(df.diamonds_mutated)



##################################################
# ARRANGE
# - arrange() sorts your data by a given variable
#
##################################################

# inspect 
head(diamonds)
min(diamonds$price) #326
max(diamonds$price) #18823


# 1: sort the data by "price" using arrange()
#    in ascending order (i.e., low to high)
#    note: arrange() sorts in descending order by default

df.diamonds_sorted <- arrange(diamonds, price)


# inspect
# note: top row has price 326
head(df.diamonds_sorted) 



# 2: sort the data by "cut" using arrange()
#    in *descending* order(i.e., high to low)
df.diamonds_sorted_desc <- arrange(diamonds, desc(price))


# inspect
# note: top row has price 18823
head(df.diamonds_sorted_desc)



#############################################################
# SUMMARIZE
# - summarize() summarizes a variable in your data frame
#   by applying a given summary statistic 
# - i.e., you provide the statistic function, and summarize
#   returns a that statistic on a particular variable
#
#############################################################

#-------------------------------
# Calculate diamond price stats
#-------------------------------

# inspect
head(diamonds)

# mean
summarize(diamonds, mean(price)) # 3932.8

# median
summarize(diamonds, median(price)) # 2401

# max
summarize(diamonds, max(price)) # 18823

# min
summarize(diamonds, min(price)) # 326

# standard deviation
summarize(diamonds, sd(price))

#--------------------------------------------------------
# NOTE: the results of these applications of summarize()
#       return the same result as the base function call.
#       For example, summarize(diamonds, mean(price))
#       returns the same value as mean(diamonds$price)
#--------------------------------------------------------




#-------------------------------------------------------------
# Calculate statistic and give the summarized variable a name
#-------------------------------------------------------------

# mean
summarize(diamonds, mean_price = mean(price)) # 3932.8

# median
summarize(diamonds, max_price = max(price)) # 2401





# Calculate maximum price-per-carat
# NOTE: this is a two step process
#       1. we first mutate the data frame to create the variable price_per_carat
#       2. we then summarize that data frame with max(price_per_carat)
df.diamonds_mutated <- mutate(diamonds, price_per_carat = price / carat)
summarize(df.diamonds_mutated, max(price_per_carat))




###################################################
# GROUP BY
# - group_by() groups your data (i.e., aggregates)
#   by a given variable
#
###################################################


#-------------------------------------------
# aggregate (i.e., group) the diamonds data
#  by the "cut" variable using group_by()
#-------------------------------------------


groupdf.diamonds_by_cut <- group_by(diamonds, cut)
class(groupdf.diamonds_by_cut)


summarise(groupdf.diamonds_by_cut, mean(price))


groupdf.diamonds_by_cut_color <- group_by(diamonds, cut, color)

summarize(groupdf.diamonds_by_cut_color, max(price))



#############################################
# min_rank()
#  - min_rank() is a helper function.
#    dplyr has several helper functions
#    and we're going to learn this one now
#    because we'll be using it later
#############################################

#--------------------------------------------
# Find 10 cheapest diamonds
#  notes: 1. to do this, we need to use 
#            a "helper function", min_rank()
#--------------------------------------------

filter(diamonds, min_rank(price) <= 10)


#-------------------------------------------------------
# Find 10 most expensive diamonds
#  notes: 1. to do this, we need to use 
#            a "helper function", min_rank()
#            AND the desc() helper function
#         2. the resulting output is sorted high-to-low
#            not low-to-high
#-------------------------------------------------------

filter(diamonds, min_rank(desc(price)) <= 10)



