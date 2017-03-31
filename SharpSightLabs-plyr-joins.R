######################################################
# JOINING DATASETS TOGETHER
#  (i.e., dplyr joins)
#
# How to join two datasets together using ...
#  - inner joins
#  - left joins
#  
#
# (Release: 2016-08)
#
# sharpsightlabs.com
# © Copyright 2016 Sharp Sight Labs
# All rights reserved
#
######################################################


#--------------
# load packages
#--------------

#library(dplyr)
library(startingDataScience)



#=============================================================
# INSPECT DATASETS
#
# In the following examples, we're going to use two datasets
#  These datasets have names of people.
#  One dataset has "first names" only. 
#  The other dataset has "first names" and "last names"
#  We will join them together on by "first name"
#=============================================================


#--------------------------
# INSPECT "FIRST NAME" data
#--------------------------

head(person_data)
str(person_data)


#----------------------------
# INSPECT "FAMILY NAME" data
#----------------------------

head(family_name_data)
str(family_name_data)



#############
# INNER JOIN
#############

inner_join(person_data, family_name_data, by = "first_name")


# first_name sex family_name
#        ned   m       stark
#     tyrion   m   lannister
#   daenarys   f   targaryen
#        rob   m       stark




#############
# LEFT JOIN
#############

left_join(person_data, family_name_data, by = "first_name")


# first_name sex family_name
#        ned   m       stark
#     tyrion   m   lannister
#   daenarys   f   targaryen
#        rob   m       stark
#      drogo   m        <NA>

