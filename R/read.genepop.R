

read.genepop <- function(filename = NULL) {
  
  genepop <- new("ARES") 
  
  # read the file into memory
  con <- file(description = filename, open = "rb")
  lines <- readLines(con)
  close(con)
  fileLength <- length(lines)
  
  #pop locations ignoring first 'title' line.
  popLocations <- grep("POP", lines[2:fileLength], ignore.case=TRUE)
  popLocations <- popLocations + 1 
  
  #resize output
  genepop@OUT <- vector("list", length(popLocations))
  
  #get title from first line
  title <- lines[1]
  
  #get Loci Column names
  lociNames <- lines[2:(popLocations[1]-1)]
  lociNames <- gsub("\t","",lociNames)
  
  #number of Loci and individuals
  npop <- length(popLocations)
  initial <- popLocations[1] 
  end <- length(lines)
  rownumber <- end - initial - npop
  genepop@number.individuals <- rownumber
  estructuras <- rep(0, rownumber)
  genepop@number.loci <- length(unlist(strsplit( lociNames, '[[:space:]]+')))
  genepop@number.pop <- npop
  
  
  #parsing of data between pops into genepop
  counter <- 1
  
  for(i in 1:length(popLocations)) {
    beginLine <- popLocations[i] + 1
    endLine <- 0
    
    group <- vector("list", length(lociNames))
    
    if(i == length(popLocations)) {	
      endLine <- length(lines) - 1 
    }
    else
    {	
      endLine <- popLocations[i+1] - 1 
    }
    
    #parse individual line
    for(line in lines[beginLine:endLine]) {
      
      #split id & alleles apart
      location <- unlist(strsplit(line, ","))
      id <- location[1]
      alleleGroup <- location[2]
      estructuras[counter] <- i
      counter <- counter + 1 
      
      #split alleles apart on whitespace
      alleles <- unlist(strsplit( alleleGroup, '[[:space:]]+'))
      alleles <- alleles[2:length(alleles)]
      
      #cat(" ",length(alleles)," => ", alleles ,"\n")	
      
      if(length(lociNames) == length(alleles)){
        
        for(j in 1:length(alleles)) 
        {
          
          #only recognize data that isn't missing
          if((alleles[j] != "00" ) && (alleles[j] != "0000") && (alleles[j] != "000000")) {
            digits <- unlist(strsplit( alleles[j], ""))
            
            if(length(digits) == 2) {
              digits1 <- unlist(substring( alleles[j], 1))
              digits2 <- unlist(substring( alleles[j], 2))
              if(digits1 != "0")
              {
                group[[j]] <- c(group[[j]], digits1)
              }
              if( digits2 != "0" )
              {
                group[[j]] <- c(group[[j]], digits2)
              }
            
            } else if(length(digits) == 4) {
              digits1 <- unlist(substring( alleles[j], 1, 2 ))
              digits2 <- unlist(substring( alleles[j], 3, 4 ))
              if( digits1 != "00" )
              {
                group[[j]] <- c( group[[j]], digits1 )
              }
              if( digits2 != "00" )
              {
                group[[j]] <- c( group[[j]], digits2 )
              }
              
            } else if(length(digits) == 6) {
            
              digits1 <- unlist(substring( alleles[j], 1, 3 ))
              digits2 <- unlist(substring( alleles[j], 4, 6 ))
              if(digits1 != "000")
              {
                group[[j]] <- c(group[[j]], digits1)
              }
              if(digits2 != "000")
              {
                group[[j]] <- c(group[[j]], digits2)
              }
            }
            
            group[[j]] <- unique(group[[j]])
            group[[j]] <- sort(group[[j]])
            
            allAttributes <- list(name=lociNames[j], id=id)
            attributes(group[[j]]) <- allAttributes
            
            
          }
          
        }
        
      }
    }
    
    attributes(group) <- list(id=i, beginLine=beginLine, endLine=endLine)
    genepop@OUT[[i]] <- group

  }
  
  contents <- list(title = title, filename = filename, lines = lines)
  
  genepop@title <- title
  genepop@lines <- lines
  
  #parsing of data between pops into genepop
  for(i in 1:length(genepop@OUT)) {
    
    group <- genepop@OUT[[i]]
    id <- attr(group,"id")
    beginLine <- attr(group,"beginLine")
    endLine <- attr(group,"endLine")
    
    colCount <- endLine - beginLine + 1
    col <- 0
    rowCount <- 0
    rowNames = list
    
    for(j in 1:length(genepop@OUT[[i]])) {
      locus <- genepop@OUT[[i]][[j]]
      rowCount <- rowCount + length(locus)
      locusName <- attr(locus,"name")
      for( k in 1:length(locus)) {
        name <- paste( list(locusName,"-",locus[k]), collapse="") 
        rowNames <- c( rowNames, name)
      }
    }
    
    rowNames <- rowNames[2:length(rowNames)]
    output <- matrix(0,nrow=rowCount, ncol=colCount, byrow=TRUE )
    
    #parse individual line
    for(line in lines[beginLine:endLine]) {
      
      #split id & alleles apart
      location <- unlist(strsplit( line, ","))
      id <- location[1]
      alleleGroup <- location[2]
      col <- col + 1
      
      #split alleles apart on whitespace
      alleles <- unlist(strsplit( alleleGroup, '[[:space:]]+'))
      alleles <- alleles[2:length(alleles)]
      
      #cat(" col( ",col," )\n")
      for( j in 1:length(alleles)) {
        #only recognize data that isn't missing
        if((alleles[j] != "0000") && (alleles[j] != "000000"))
        {
          digits <- unlist(strsplit( alleles[j], ""))
          
          if(length(digits) == 4) {
            digits1 <- unlist(substring( alleles[j], 1, 2 ))
            digits2 <- unlist(substring( alleles[j], 3, 4 ))
          }
          else if( length(digits) == 6 )
          {
            digits1 <- unlist(substring(alleles[j], 1, 3))
            digits2 <- unlist(substring(alleles[j], 4, 6))
          }
          
          existDigits1 <- grep(digits1, group[[j]])
          existDigits2 <- grep(digits2, group[[j]])
          
          if(length(existDigits1) > 0) {
            matchName <- paste( list(attr(group[[j]],"name"),"-",digits1), collapse="") 
            rowNum <- grep(matchName, rowNames)
            if( length(rowNum) > 0 ) {
              output[rowNum,col] <- output[rowNum,col] + 1
            }
          }
          
          if(length(existDigits2) > 0) {
            matchName <- paste( list(attr(group[[j]],"name"),"-",digits2), collapse="") 
            rowNum <- grep(matchName, rowNames)
            if( length(rowNum) > 0) {
              output[rowNum,col] <- output[rowNum,col] + 1
            } 
          }
        }
      }
      
      attr(genepop@OUT[[i]], "output") <- output
      
    }
  }
  
  # output
  genepop@structures <- factor(estructuras)
  genepop
}
