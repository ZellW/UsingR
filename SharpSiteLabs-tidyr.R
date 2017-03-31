#############################################################
# DATA RESHAPING WITH TIDYR
#
# How to ...
#   - reshape data wide to long:    spread()
#   - reshape data long to wide:    gather()
#   - split a variable:             separate()
#   - combine two varibale values:  unite()
#
# Also: more examples at the bottom
#       - specifically, examples relevant to machine learning
#
# (Release: 2016-08)
#
# sharpsightlabs.com
# © Copyright 2016 Sharp Sight Labs
# All rights reserved
#
#############################################################



#----------------
# LOAD LIBRARIES
#----------------

#library(tidyr)
#library(tempPackage)
library(startingDataScience)


#################
# tidyr functions
#################



#=================================
# GATHER
# - gather() converts wide to long
#
#=================================

# inspect
head(revenue_2_untidy)


# gather()
# - reshape wide to long
# - new key: year
# - new value: income

gather(revenue_2_untidy, quarter, revenue, -region)




#===================================
# SPREAD
# - spread() converts long to wide
#
#===================================

# inspect
print(revenue_tidy)


# spread()
# - reshape long to wide
# - "year" variable holds the values that will become new columns

spread(revenue_tidy , quarter, revenue)




#====================================================
# SEPARATE
# - separate() divides a column into multiple columns
#
#====================================================

# inspect
print(founders)


# separate()
# - break "name" into "first_name" and "last_name"
separate(founders , name , c("first_name","last_name") , sep = " ")




#====================================================
# UNITE
# - unite() combines multiple columns into one column
#
#====================================================


# inspect
print(founders_2)


# unite()
# - combine
unite(founders_2, name, first_name, last_name, sep = " ")






################
# MORE EXAMPLES
################


#-----------------------------------------------------
# EXAMPLE:
#  - Boston dataset
#  - 'Boston' is a good machine learning dataset
#    and what we're doing here is typical ML
#    data exploration
#  - we're going to create a density plot of each
#    variable (all in the same plot)
#  - to do this, we'll need to reshape
#-----------------------------------------------------


# get data
data(Boston, package = "MASS")


# inspect
head(Boston) # print first few rows
ncol(Boston) # calculate number of columns


# gather()
# - reshape wide to long
# - we want to have our variable names along the rows

df.boston_gathered <- gather(Boston)

# inspect
head(df.boston_gathered)


# plot density plot of each variable
ggplot(data = df.boston_gathered, aes(x = value)) +
  geom_density() +
  facet_wrap(~ key, scales = "free")


# REDO using dplyr pipes:
# - here, we'll redo the process of reshaping and plotting Boston
#   ... we'll use dplyr pipes to do it all-in-one

Boston %>%
  gather() %>%
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~ key, scales = "free")




#----------------------------------------------------------------------
# EXAMPLE:
# - bbbDescr dataset from 'caret'
# - note: again, this is a machine learning dataset and
#         what we're doing here is typical of ML pre-modeling work
# - ultimately, we want to create density plots of all of our variables
# - we need to reshape to get the data in a form that ggplot can use
#   in a small-multiple chart (i.e., faceting)
#----------------------------------------------------------------------


# get data from 'caret' package
library(caret)
data(BloodBrain)


# inspect
head(bbbDescr)  # print top few rows
dim(bbbDescr)   # get number of rows & columns ... that's a lot!




# gather()
# - reshape wide to long
# - we want to have our variable names along the rows

df.bbb_gathered <- gather(bbbDescr)


# inspect
# - note: all our variables are allong the rows now
head(df.bbb_gathered)


# plot small multiple
# density plot of each separate variable

ggplot(data = df.bbb_gathered, aes(x = value)) +
  geom_density() +
  facet_wrap(~ key, scale = "free")






