######################################################
# The dplyr Pipe Operator
#  AKA, "chaining"
#
# the '%>%' operator allows you to "chain" together
# multiple commands, function, and dplyr verbs
# (this is somtimes called 'piping', much like pipes
#  in UNIX)
#
#
# (Release: 2016-08)
#
# sharpsightlabs.com
# © Copyright 2016 Sharp Sight Labs
# All rights reserved
#
######################################################



#---------------
#LOAD LIBRARIES
#---------------

library(dplyr)
library(ggplot2)



#------------------
# DATA INSPECTION
#------------------

# get first few records
head(diamonds)

# examine the levels of the "cut" variable
levels(diamonds$cut)

# create a frequency table of the "cut" variable
table(diamonds$cut)

# create a *proportion* table of the "cut" variable
table(diamonds$cut) %>%
  prop.table()




# CREATE DATA WITH ONLY 'carat' & 'price' variables
diamonds %>%
  select(carat, price)



###########
# FILTER
###########

#-----------------------------
# Filter
# - 'Ideal' cut diamonds only
#-----------------------------

diamonds %>%
  filter(cut == 'Ideal')


# maximum price for 'Ideal' cut diamonds
diamonds %>%
  filter(cut == 'Ideal') %>%
  summarize(max(price))


# create subset of 'Fair' cut diamonds
diamonds %>%
  filter(cut == 'Fair')


# Create histogram of price 
#  for only 'Fair' cut diamonds
diamonds %>%
  filter(cut == 'Fair') %>%
  ggplot(aes(x = price)) +
    geom_histogram()


# Create boxplot of carat, by color variable
#  for subset of 'Fair' cut diamonds
diamonds %>%
  filter(cut == 'Fair') %>%
  ggplot(aes(x = color, y = carat)) +
    geom_boxplot()


# REDO: 
#  - Create boxplot of carat, by color variable
#   for subset of 'Fair' cut diamonds
#  - Flip x/y coordinates
diamonds %>%
  filter(cut == 'Fair') %>%
  ggplot(aes(x = color, y = carat)) +
    geom_boxplot() +
    coord_flip()



# Calculate maximum price for 'E' color diamonds
#  grouped by different diamond cuts
diamonds %>%
  filter(color == 'E') %>%
  group_by(cut) %>%
  summarize(max_price = max(price))



# Make a bar chart of maximum price for 'E' color diamonds
#  grouped by different diamond cuts
diamonds %>%
  filter(color == 'E') %>%
  group_by(cut) %>%
  summarize(max_price = max(price)) %>%
  ggplot(aes(x = cut, y = max_price)) +
    geom_bar(stat = "identity")



# REDO:
# - Make a bar chart of maximum price for 'E' color diamonds
#   grouped by different diamond cuts
#   AND flip the coordinates to create a horizontal bar chart

diamonds %>%
  filter(color == 'E') %>%
  group_by(cut) %>%
  summarize(max_price = max(price)) %>%
  ggplot(aes(x = cut, y = max_price)) +
    geom_bar(stat = "identity") +
    coord_flip()





#############
# SELECT
#############

diamonds %>%
  select(carat, cut, price)
  


#--------------------------------------
# Reshape 'diamonds' dataset
# - reduce variables
# - add new varibale 'price_per_carat'
#--------------------------------------


# select carat, cut, and price
#  then create a new variable 'price_per_carat'
diamonds %>%
  select(carat, cut, price) %>%
  mutate(price_per_carat = price / carat)



# Redo: pipe into head
diamonds %>%
  select(carat, cut, price) %>%
  mutate(price_per_carat = price / carat) %>%
  head()



# Calculate maximum price per carat
#  for different diamond cuts
diamonds %>%
  select(carat, cut, price) %>%
  mutate(price_per_carat = price / carat) %>%
  group_by(cut) %>%
  summarize(max(price_per_carat))



diamonds %>%
  select(carat, cut, price) %>%
  mutate(price_per_carat = price / carat) %>%
  group_by(cut) %>%
  summarize(max_price_per_carat = max(price_per_carat)) %>%
  ggplot(aes(x = cut, y = max_price_per_carat)) +
    geom_bar(stat = "identity")



diamonds %>%
  select(carat, cut, price) %>%
  mutate(price_per_carat = price / carat) %>%
  group_by(cut) %>%
  summarize(max_price_per_carat = max(price_per_carat)) %>%
  ggplot(aes(x = cut, y = max_price_per_carat)) +
    geom_bar(stat = "identity") +
    coord_flip()











######################
# SUMMARIZE by GROUPS
######################


#------------------------------
# Average diamond price by cut
#  (diamonds example)
#------------------------------

# EXAMPLE WITHOUT CHAINING
group.diamonds_by_cut <- group_by(diamonds, cut)
summarise(group.diamonds_by_cut, avg_price = mean(price))  


# REDO _WITH_ CHAINING
diamonds %>%
  group_by(cut) %>%
  summarise(mean(price))


#      Fair    4358.758
#      Good    3928.864
# Very Good    3981.760
#   Premium    4584.258
#     Ideal    3457.542


#--------------------------------
# Median diamond weight (carats)
#  by cut
#--------------------------------
diamonds %>%
  group_by(cut) %>%
  summarise(median(carat))





############################
# MORE COMPLICATED EXAMPLES
############################







# CREATE DATA WITH ONLY 'carat' & 'price' variables
# AND PIPE INTO ggplot()

diamonds %>%
  select(carat, price) %>%
  ggplot(aes(x = carat, y = price)) +
    geom_point()





diamonds %>%
  filter(cut == 'Ideal') %>%
  ggplot(aes(x = carat, y = price)) +
    geom_point()


  

#---------
# ARRANGE
#---------

# sort by price, low to high
diamonds %>%
  arrange(price)


# find the cheapest diamonds
diamonds %>%
  filter(min_rank(price) >= 10)
  

# find the cheapest diamonds
diamonds %>%
  filter(min_rank(price) >= 10) %>%
  arrange()



# sort by price, high to low
diamonds %>%
  arrange(desc(price))




  