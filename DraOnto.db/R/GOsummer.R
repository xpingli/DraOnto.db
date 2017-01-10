#' GOsummer is the function that you can use to do summerization on the myGO object.
#' 
#' @param mygo,Type,Term
#' 
#' Type is the ontology type: MF, BP, CC, PC
#' Term is the GOterms or terms you are interested to see, such as response, binding, et.
#' 
#' MF terms include: "transporter activity", "binding", "catalytic activity", "structural molecular activity", "receptor activity", "antioxidant activity"
#' 
#' BP terms include:"cellular process", "localization", "metabolic process", "response to stimulus", "biological regulation", "cellular process", "cellular component organization or biogenesis"
#' 
#' CC terms include: "cell part", "marcromolecular complex", "membrane", "organelle"
#' 
#' PC terms include:"chaperone", "hydrolase","transporter", "oxidoreductase", "ligase", "enzyme modulator", "isomerase", "nucleic acid binding", "receptor", "transcription factor", "transfer/carrier protein", "transferase", "defense/immunity protein", "ribosomal RNA",  "tRNA"
#' 
#' IF inputs are not in Type or Term, it returns a list of terms you can choose from
#' 
#' IF Type input is not in the type vector, and you know the Term, it returns a list of summaries based on the Term
#' 
#' IF Term is unknown, and you know the Type, it returns a list of summaries based on the Type
#' 
#' IF both inputs are in the function, it returns a list of summaries according to the Type and Term
#' 
#' @examples 
#' 
#' \dontrun{GOsummer(mygo.e.34hr, "MF", "transporter activity")}
#' \dontrun{GOsummer(mygo.e.34hr, "I don't know", "transporter")}
#' @export
GOsummer <- function(mygo, Type, Term){
        
        types <- c("MF", "BP", "CC", "PC")
        
        mfterms <- c("transporter activity", "binding", "catalytic activity", "structural molecular activity", "receptor activity", "antioxidant activity")
        bpterms <- c("cellular process", "localization", "metabolic process", "response to stimulus", "biological regulation", "cellular process", "cellular component organization or biogenesis")
        ccterms <- c("cell part", "marcromolecular complex", "membrane", "organelle")
        pcterms <- c("chaperone", "hydrolase","transporter", "oxidoreductase", "ligase", "enzyme modulator", "isomerase", "nucleic acid binding", "receptor", "transcription factor", "transfer/carrier protein", "transferase", "defense/immunity protein", "ribosomal RNA",  "tRNA")
        
        terms <- unique(unlist(strsplit(unique(as.character(mygo$goTerm)), " ")))
        
        
        
        if(! Type %in% types & ! Term %in% terms){
                
                warning(cat("Types are:", types, sep = ","))
                warning(cat("\nMF goterms:", mfterms, sep = ","))
                warning(cat("\nBP goterms:", bpterms, sep = ","))
                warning(cat("\nCC goterms:", ccterms, sep = ","))
                warning(cat("\nPC goterms:", pcterms, sep = ","))
                
        } else if( Type %in% types & ! Term %in% terms ){
                
                filtered <- mygo[mygo$type == Type, ]
                
                filtered <- droplevels(filtered)
                
                # check up and down regulation
                
                table <- table(filtered$goTerm, regulation = ifelse(filtered$DE == -1, "down", "up"))
                
                
                upgenes <- filtered[filtered$DE == 1,]
                downgenes <- filtered[filtered$DE == -1,]
                
                ups <- sum(filtered$DE == 1)
                downs <- sum(filtered$DE == -1)
                
                upIDs <- split(upgenes$refID, upgenes$goTerm)
                downIDs <- split(downgenes$refID, downgenes$goTerm)
                
                
                ## output
                print(filtered)
                cat("\n")
                print(paste("up regulated genes:", ups))
                print(paste("down regulated genes:", downs))
                cat("\n")
                print(table)
                cat("\n")
                print("Genes upreguated:")
                print(upIDs)
                cat("\n")
                print("Genes downreguated:")
                print(downIDs)
                
                
        } else if (! Type %in% types & Term %in% terms){
                
                patt_index <- grep(paste("*", Term, "*", sep = ""), mygo$goTerm, value = FALSE)
                
                filtered <- mygo[patt_index, ]
                
                filtered <- droplevels(filtered)
                
                df_split <- split(filtered, filtered$type)
                
                nam <- vector()
                ups <- vector()
                upids <- list()
                downs <- vector()
                downids <- list()
                biglist <- list()
                for(i in 1:length(df_split)){
                        
                        nam[i] <- names(df_split)[i]
                        TT <- df_split[[i]]
                        ups[i] <- sum(TT$DE == 1)
                        downs[i] <- sum(TT$DE == -1)
                        
                        upids[[i]] <- as.character(TT[TT$DE == 1, ]$refID)
                        downids[[i]] <- as.character(TT[TT$DE == -1,]$refID)
                        biglist[[i]] <- list(goType = nam[i], up_refID = upids[[i]], down_refID = downids[[i]])
                        
                        
                        
                }
                
                df_term <- data.frame(GOtype = nam, Up = ups, Down = downs)
                
                print(df_term)
                cat("\n")
                print(biglist)
                
                
        } else if( Type %in% types & Term %in% terms) {
                go_df <- subset(mygo , type == Type & mygo$goTerm == Term)
                
                go_df <- droplevels(go_df)
                
                table <- table(go_df$goTerm, ifelse(go_df$DE == -1, "down", "up"))
                
                upgenes <- go_df[go_df$DE == 1, ]
                downgenes <- go_df[go_df$DE == -1,]
                
                
                upID <- upgenes$refID
                downID <- downgenes$refID
                
                
                #output
                print(go_df)
                cat("\n")
                print(table)
                cat("\n")
                print(paste("Up regulated genes:", length(upID))) 
                cat("\n")
                print(paste( "Down regulated genes:", length(downID)))
                cat("\n")
                print(paste("up IDs:", upID))
                cat("\n")
                print(paste("down IDs:", downID))
        }
        
}
