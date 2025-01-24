powo_distribution <- function(sp){
  sp <- gsub("_"," ",sp)
  powo <- search_powo(sp)
  powo_info <- tidy(powo)
  
  distribution <- NULL
  no_match <- F
  
  require(plyr)
  
  if (nrow(powo_info) == 0) {
    message("no match found")
    no_match <- T
    tdwg3_code <- data.frame(sp=sp,tdwg3_code=NA)
    return(tdwg3_code)
  }
  
  accepted_match <- subset(powo_info,(name == paste0(strsplit(sp," ")[[1]][1], " × ",strsplit(sp," ")[[1]][2]) | name == sp) & accepted == T) #x for Erysimum cheiri??
  
  if (nrow(accepted_match) == 0) {
    message("no match found")
    no_match <- T
    tdwg3_code <- data.frame(sp=sp,tdwg3_code=NA)
    return(tdwg3_code)
  }
  
  synonym <- 1
  
  if (!("synonymOf" %in% colnames(powo_info))) {
    synonym <- 0
    message("no synonym found")
    powo_info <- subset(powo_info,(name == paste0(strsplit(sp," ")[[1]][1], " × ",strsplit(sp," ")[[1]][2]) | name == sp) & accepted == T)
    powo_ID <- sub(".*names:","",powo_info$fqId)
  }
  
  if ("synonymOf" %in% colnames(powo_info)) {
    if (nrow(accepted_match) != 0 & is.null(nrow(accepted_match$synonymOf))) {
      synonym <- 0
      message("found synonym, but the name is correct")
      powo_ID <- sub(".*names:","",accepted_match$fqId)
    }
  }
  
  if ("synonymOf" %in% colnames(powo_info)) {
    if (!is.null(powo_info[powo_info$name == sp,"synonymOf"]) & synonym == 1) {
      message("found synonym")
      synonym <- do.call(rbind.fill,powo_info$synonymOf)
      powo <- lapply(synonym$name,search_powo)
      powo_info <- do.call(rbind.fill,lapply(powo,tidy))
      powo_info <- powo_info[powo_info$name %in% synonym$name & powo_info$accepted == T,]
      powo_ID <- sub(".*names:","",powo_info$fqId)
    } 
  }
  
  tdwg3_code_all <- NULL
  for (j in 1:length(powo_ID)){
    record <- lookup_powo(powo_ID[[j]], distribution=TRUE)
    tidied <- tidy(record)
    
    if (!"natives" %in% colnames(as.data.frame(tidied$distribution))){
      message("no native distribution")
      tdwg3_code <- NA
      next()
    }  
    
    if ("distribution" %in% colnames(tidied)) {
      tidy_distribution <- tidied %>%
        select(fqId, distribution) %>%
        unnest(cols=distribution) %>%
        unnest(cols=natives)
      
      tdwg3_code <- tidy_distribution$tdwgCode
    }
    tdwg3_code_all <- rbind(tdwg3_code_all,data.frame(sp=sp,tdwg3_code=tdwg3_code))
  }

  return(tdwg3_code_all)
}
