#####################################################
# LUBRIDATE
# Manipulating dates and times
#
# How to ...
#   - parse character data into dates
#   - access date attributes
#   - modify/change date attributes
#
#
# (Release: 2016-11)
#
# sharpsightlabs.com
# © Copyright 2016 Sharp Sight Labs
# All rights reserved
#
#####################################################



library(startingDataScience)



###############################################
# DATE PARSING FUNCTIONS
# 
# ymd()
# ydm()
# dmy()
# dym()
# mdy()
# myd()
# 
# NOTE: - these functions can parse dates that 
#         are in a variety of structures
#       - note that in the following examples, 
#         the "separators" are different 
#
###############################################


#---------------
# YEAR MONTH DAY
#  ymd()
#---------------

# parse date
date.bf <- ymd("1706-01-17")


# inspect value
print(date.bf)


# inspect class
class(date.bf)

#---------------
# DAY MONTH YEAR
#  dmy()
#---------------

# parse date
date.nt <- dmy("10/07/1856")


# inspect value
print(date.nt)


# inspect class
class(date.nt)



#-----------------------
# PARSING MULTIPLE DATES
#-----------------------

# store vector of two strings 
vec.two_strings <- c("1856/07/10","1706/01/17")


# parse the strings into proper dates
vec.two_dates <- ymd(vec.two_strings)


# inspect
print(vec.two_dates)




################################################################
# ACCESSING DATE ATTRIBUTES
#  To access individual date attribues (i.e., day, month, year)
#  we can use a set of "accessor funcitons
# 
# year()    : year access
# month()   : month
# week()    : week
# wday()    : weekday
# day()     : day
# hour()    : hour
# minute()  : minute
# second()  : second 
# tz()      : time zone
#
################################################################


#-----------------------
# ACCESS DATE ATTRIBUTES
#-----------------------

# parse date
date.charlie <- mdy("January 1, 1924")


# inspect
date.charlie


# access individual date attributes
year(date.charlie)
month(date.charlie)
day(date.charlie)


#-----------------------
# access time attributes
#-----------------------


# get current time
#  now()
datetime.now <- now("Singapore")


#inspect
print(datetime.now)


# retrieve time attributes
hour(datetime.now)
minute(datetime.now)
second(datetime.now)




###########################################
# CHANGING DATE ATTRIBUTES
#  To change date attributes, we can use 
#  the lubridate "accessor functions"
# 
###########################################


#------------------------
# SETTING DATE ATTRIBUTES
#------------------------

# create "dummy" date
date.wb <- ymd("1900-01-01")


# inspect (retrieve individual attributes)
year(date.wb)
month(date.wb)
day(date.wb)


# assign new values for date attributes
year(date.wb) <- 1930
month(date.wb) <- 8
day(date.wb) <- 30


# inspect again
#  note that the date is different
print(date.wb)



#--------------------------------
# SET MULTIPLE ATTRIBUTES
# update(): setting multiple date 
#           attributes at once
#--------------------------------

# set a "dummy" date
date.rf <- ymd("1900-01-01")


# "update" year, month, and day to correct values
date.rf <- update(date.rf, year = 1918, month = 5, day = 11)


# inspect
print(date.rf)



#------------------------
# SETTING TIME ATTRIBUTES
#------------------------

# create date-time variable
datetime.now <- now("Singapore")


# inspect
print(datetime.now)


# assign new values for time attributes
hour(datetime.now) <- 0
minute(datetime.now) <- 0
second(datetime.now) <- 42


# re-inspect
print(datetime.now)


