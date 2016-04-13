
setClass("ARES", 
         
         representation(OUT = "list",
                        structures = "factor",
                        number.pop = "numeric",
                        number.individuals = "numeric", 
                        number.loci = "numeric", 
                        title = "character",
                        lines = "character")
        )

setMethod("show", "ARES", function(object) {
  cat("\n", object@title, "\n")
  npop <- object@number.pop
  cat(" Number of populations:", npop, "populations", "\n")
  cat(" Populations codes:", 1:npop, "\n")
  cat(" Number of loci:", object@number.loci, "\n")
  cat(" Number of individuals:", object@number.individuals, "\n\n")
  estructuras <- object@structures
  ones <- rep(1, length(estructuras))
  indpp <- tapply(ones, estructuras, sum)
  subpop <- data.frame(1:object@number.pop, indpp)
  colnames(subpop) <- c("pop", "IndNum")
  cat(" Structures:", "\n") 
  print(subpop) 
  cat("\n")
})
