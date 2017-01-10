#' A myGO function
#' 
#' This allows you to do Deinococcus radiodurans ontology analysis. The funciton generates a data frame with gene symbol, refID, gene ID, Differential expression pattern, annotation from NCBI, GO terms from Panther.
#' 
#' The data base is based on all the genes that have hits in PATHNER from our previous analysis
#' 
#' You need MF, BP, CC, PC for database construction for this function to run
#' 
#' You need to input the annotation dataset from annotation.rmd
#' 
#' write = T to write a csv file with all the genes in the ontology database
#' 
#' @examples
#' \dontrun{myGO(MF, BP, CC, PC, edgeR_34hr, FALSE)}
#' 
#' @export
myGO <- function(MF, BP, CC, PC, annoset, write){
        mf <- MF[,-1]
        bp <- BP[,-1]
        cc <- CC[,-1]
        pc <- PC[,-1]
        #colname consistence
        colN <- colnames(pc)
        colnames(mf) <- colN
        colnames(bp) <- colN
        
        #remove duplicate rows
        
        d_mf <- duplicated(mf)
        d_bp <- duplicated(bp)
        d_cc <- duplicated(cc)
        d_pc <- duplicated(pc)
        
        
        
        mf2 <- mf[!d_mf,]
        bp2 <- bp[!d_bp,]
        cc2 <- cc[!d_cc,]
        pc2 <- pc[!d_pc,]
        
        # global ontology dataset
        combined2 <- rbind(mf2, bp2, cc2, pc2)
        duplicat <- duplicated(combined2)
        
        
        # cleaned ontology dataset
        GO.dra.db <- combined2[!duplicat,]
        
        if(write == TRUE){
                write.csv(combined2, "[dra]withduplicate.goterm.csv")
                write.csv(GO.dra.db, "[GO]dra.database.csv")
        } else {
                warning("Not to produce the ontology library in the same filepath.")
        }
        
        
        reGO<- function(annoset, goterm){
                
                anno <- read.csv(annoset, stringsAsFactors = FALSE, strip.white = T)[,-1]
                anno <- as.data.frame(apply(anno, 2, trimws, "both"), stringsAsFactors = FALSE)
                anno.refID <- anno$refseq_locus
                
                go <- as.data.frame(apply(goterm, 2, as.character), stringsAsFactors = FALSE)
                
                go <- as.data.frame(apply(go, 2, trimws, "both"), stringsAsFactors = FALSE)
                
                go[is.na(go$symbol),]$symbol <- "No Match"
                
                go.refID <- go$refID
                
                # from anno set
                geneid <- vector()
                dexpr <- vector()
                notation <- vector()
                
                # match first
                k <- go.refID %in% anno.refID
                set <- go[which(k == TRUE), ]
                
                for (i in 1:nrow(set)){
                        
                        #from annotation dataset
                        geneid[i] <- anno[anno$refseq_locus == set$refID[i],]$geneID
                        dexpr[i] <- anno[anno$refseq_locus == set$refID[i],]$de
                        notation[i] <- anno[anno$refseq_locus == set$refID[i],]$annotation
                }
                
                
                df <- as.data.frame(cbind(  symbol = set$symbol,
                                            refID = set$refID, 
                                            geneID = geneid,
                                            DE = dexpr,
                                            annotation = notation,
                                            type = set$type,
                                            goTerm = set$Term))
                #clean the duplicates
                dup <- duplicated(df)
                df <- df[!dup,]
                df
                
        }
        
        goterm <- GO.dra.db
        
        mygo <- reGO(annoset, goterm)
        
        mygo
}






