library('RSQLite')
if(!file.exists("/home/ratanond/Desktop/Masters_Project/CCA/Tool/tool.db")){
  stop("Database does not exist!")
}

db = dbConnect(SQLite(), dbname = "/home/ratanond/Desktop/Masters_Project/CCA/Tool/tool.db")
writeLines("Starting Sparse Canonical Correlation Analysis")

while(TRUE) {
res <- dbGetQuery(db, "SELECT * from ccajobs where (status=='incomplete')")
if(nrow(res) > 0){

    row <- dbGetQuery(db, "SELECT name, project, rand, id from ccajobs where (status=='incomplete') limit 1")
    #args = c(row[1,1],row[1,2])
    randname = paste(row[1,3],row[1,4], sep="")
    writeLines(paste("Processing job",randname,sep=": "))


result <- tryCatch({
# Use print.eval=T to get plots output  
source("/home/ratanond/Desktop/Masters_Project/CCA/Tool/sCCA.R", echo = T, print.eval = T)
write.csv(x = ccaScores_old, file = outfile)

dbSendQuery(conn = db, sprintf("update ccajobs set status='Complete' where rand='%s'",row[1,3]))
}, error = function(err) {
	dbSendQuery(conn = db, sprintf("update ccajobs set status='Errored' where rand='%s'",row[1,3]))

})

## close bracket for if(res)
}
Sys.sleep(1)
## close bracket for while(TRUE) loop
}

dbDisconnect(db)
