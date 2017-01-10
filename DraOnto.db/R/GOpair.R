#' GOpair is a function to compare two mygo objects. it compares the common genes and also returns unique genes for the object
#' 
#' @export
GOpair <- function(mygo1, mygo2){
        
        
        hits1 <- length(unique(mygo1$refID))
        hits2 <- length(unique(mygo2$refID))
        
        comm1 <- mygo1$refID %in% mygo2$refID
        comm2 <- mygo2$refID %in% mygo1$refID
        
        dat1 <- mygo1[which(comm1 == TRUE),c(2,4,6,7)]
        dat1 <- droplevels(dat1)
        table1 <- table(dat1$goTerm, regulation = ifelse(dat1$DE == -1, "down", "up"))
        
        dat2 <- mygo2[which(comm2 == TRUE), c(2,4,6,7)]
        dat2 <- droplevels(dat2)
        table2 <- table(dat2$goTerm, regulation = ifelse(dat2$DE == -1, "down", "up"))
        
        uniq1 <- mygo1[which(comm1 == FALSE), c(2,4, 6,7)]
        uniq1 <- droplevels(uniq1)
        table3 <- table(uniq1$goTerm, regulation = ifelse(uniq1$DE == -1, "down", "up"))
        
        uniq2 <- mygo2[which(comm2 == FALSE), c(2,4,6,7)]
        uniq2 <- droplevels(uniq2)
        table4 <- table(uniq2$goTerm, regulation = ifelse(uniq2$DE == -1, "down", "up"))
        
        similarity1 <- mean(comm1)*100
        similarity2 <- mean(comm2)*100
        
        summary1 <- cat("Ontology hits for the first dataset:", hits1, "\nOntology hits for the second dataset:", hits2, "\nfirst data set has similarity of", round(similarity1,2), "% with the second data set", "\nSecond data set has similarity of", round(similarity2, 2), "% with the first one.")
        cat("\n")
        print(summary1)
        cat("\n")
        cat("\n")
        if(similarity1 < similarity2 & similarity2 == 100){
                print("Common genes:")
                print(dat1)
                cat("\n")
                print("expression pattern in the common genes")
                print(table1)
                cat("\n")
                print("Unique genes in mygo set 1")
                print(uniq1)
                cat("\n")
                print("expression pattern in the unique genes in mygo set1")
                print(table3)
                
                
        }else if(similarity1 > similarity2 & similarity1 == 100){
                print("Common genes:")
                print(dat2)
                cat("\n")
                print("expression pattern in the common genes")
                print(table2)
                cat("\n")
                print("Unique genes in mygo set 2")
                print(uniq2)
                cat("\n")
                print("expression pattern in the unique genes in mygo set2")
                print(table4)
                
        }else if(!similarity1 == 100 & ! similarity2== 100){
                print("Common genes:")
                print(dat1)
                cat("\n")
                print("expression pattern in the common genes")
                print(table1)
                cat("\n")
                print("Unique genes in mygo set 1")
                print(uniq1)
                cat("\n")
                print("expression pattern in the unique genes in mygo set1")
                print(table3)
                cat("\n")
                print("Unique genes in mygo set 2")
                print(uniq2)
                cat("\n")
                print("expression pattern in the unique genes in mygo set2")
                print(table4)
        }
        
        
        
        
        
} 