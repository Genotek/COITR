library(optparse)
library(RJSONIO)
library(stringr)
library(e1071)
library(SDMTools)

option_list <- list( 
  make_option("--ind1", action = "store", default = NULL, type = 'character', help = "json-file of the first parent"), 
  make_option("--ind2", action = "store", default = NULL, type = 'character', help = "json-file of the second parent"), 
  make_option("--phen", action = "store", default = "Phenotypes.json", type = 'character', help = "probabilities"),
  make_option("--out", action = "store", default = "output_all/", type = 'character', help = "folder for results") 
)

args <- parse_args(OptionParser(option_list = option_list)) 

source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep = "/"))
}
source_local('functions.r')

#get probabilities
phen <- fromJSON(args$phen)

#get parents' genotypes
ind1 <- fromJSON(args$ind1)
ind2 <- fromJSON(args$ind2)


#result for child
result_ch <- vector(mode = "list", length = 3)
names(result_ch) <- c("eyes","hair","skin")

result_ch <- lapply(names(result_ch), function(x){
  #list of snps
  snps <- names(phen[[x]])[substring(names(phen[[x]]), 1, 2) %in% "rs" | substring(names(phen[[x]]), 1, 3) %in% "gen"]
  mom <- rep(NA, length(snps))
  names(mom) <- snps
  mom[names(ind1[[x]]$Snps)] <- ind1[[x]]$Snps 
  mom = sapply(mom,drop.null)
  mom[which(mom == "NA")] <- NA
  dad <- rep(NA, length(snps))
  names(dad) <- snps
  dad[names(ind2[[x]]$Snps)] <- ind2[[x]]$Snps 
  dad = sapply(dad,drop.null)
  dad[which(dad == "NA")] <- NA
  res_mid <- child_pred(mother = mom, father = dad, snps = snps, trait = x, phenes = phen)
  result_ch[[x]] <- stat_child(comb_prob = res_mid$comb_prob, comb_prob_2 = res_mid$comb_prob_2, phenes = phen, trait = x)
})

names(result_ch) <- c("eyes","hair","skin")

exportJson <- toJSON(result_ch)
write(exportJson, file = paste0(args$out,"/","child.json"))
