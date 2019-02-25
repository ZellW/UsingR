#A remote database to play with: 
###Server Name: msedxeus.database.windows.net 
###Database: DAT209x01 
###Login: RLogin 
###Password: P@ssw0rd

library(RODBC)
#Connect to database 
connStr <- paste(
  "Server=msedxeus.database.windows.net", 
  "Database=DAT209x01", "uid=Rlogin", 
  "pwd=P@ssw0rd", "Driver={SQL Server}", sep=";") 
conn <- odbcDriverConnect(connStr)

#Use sqlTables to list all tables in the database. 
#Submit query from R :
tab <- sqlTables(conn)
head(tab)

#Use sqlFetch to get a table from the database. 
#Get the table ’manufacturer’ from SCHEMA ’bi’:
mf <- sqlFetch(conn,"bi.manufacturer")
mf

#Use sqlQuery for more advanced queries. 
#SQL syntax example: 
#SELECT Manufacturer FROM bi.manufacturer WHERE ManufacturerID < 10
#Submit query to R :
query <- "SELECT Manufacturer FROM bi.manufacturer WHERE ManufacturerID < 10"
sqlQuery(conn, query)

#A common use case: Fetching entire table is infeasible
#Get some info without complete fetch:
##Count number of rows in table ’salesFact’: 
sqlQuery(conn, "SELECT COUNT(*) FROM bi.salesFact")
##Show some column info
sqlColumns(conn,"bi.salesFact")[c("COLUMN_NAME","TYPE_NAME")]
##Show the first 2 rows:
sqlQuery(conn, "SELECT TOP 2 * FROM bi.salesFact")
##Fetch a subset
df <- sqlQuery(conn, "SELECT * FROM bi.salesFact WHERE Zip='30116'")
dim(df)
##Classes of variables on the R side
sapply(df, class)
#(Recall that the variable ’Zip’ was stored as the SQL speciﬁc type ’varchar’. When read into R it became an integer.
# To avoid conversion you may want to pass as.is=TRUE to your SQL query.)

#random example
df <- sqlQuery(conn, "SELECT AVG(Revenue), STDEV(Revenue), Zip FROM bi.salesFact GROUP BY Zip")
colnames(df) <- c("AVG(Revenue)", "STDEV(Revenue)", "Zip")
head(df)

#Close the connection
close(conn)