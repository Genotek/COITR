library(optparse)
library(RJSONIO)
library(stringr)
library(e1071)
library(SDMTools)

option_list <- list(  
  make_option("--txt", action = "store", default = NULL, type = 'character', help = "chip txt file"), 
  make_option("--missed", action = "store", default = "missed/aa0000.txt.missed", type = 'character', help = "missed file"),
  make_option("--all", action = "store", default = "population.csv", type = 'character', help = "all possible"),
  make_option("--phen", action = "store", default = "Phenotypes.json", type = 'character', help = "probabilities for phenotypes"),
  make_option("--out", action = "store", default = "output_all/", type = 'character', help = "folder for results") 
)

args <- parse_args(OptionParser(option_list = option_list)) 

source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
}
source_local('functions.r')

chip <- read.table(args$txt, col.names = c('id', 'chr', 'pos', 'type'), colClasses = c("character","factor","integer","factor")) 
if (file.exists(args$missed)){
  missed <- read.table(args$missed, col.names = c('id', 'chr', 'pos', 'type', 'ch_1', 'ch_2'), colClasses = c("character","factor","integer","factor","NULL","NULL"))
  chip <- rbind(chip,missed)
}

#get probabilities
phen <- fromJSON(args$phen)

result <- vector(mode="list", length=5)
names(result) <- c("eyes","hair","skin", "frec", "child")

#determine person phenotypes
result$eyes <- det_phen(ph = phen$eyes,  chip = chip)
result$hair <- det_phen(ph = phen$hair,  chip = chip)
result$skin <- det_phen(ph = phen$skin,  chip = chip)
result$frec <- det_phen(ph = phen$frec,  chip = chip)

#the genotypes of population
population <- read.csv(args$all, header = T)[,-1]

#result for child
result_ch = lapply(names(result)[1:3], function(x){
  snps <- names(phen[[x]])[substring(names(phen[[x]]), 1, 2) %in% "rs" | substring(names(phen[[x]]), 1, 3) %in% "gen"]
  mot <- rep(NA, length(snps))
  names(mot) <- snps
  mot[names(result[[x]]$Snps)] <- result[[x]]$Snps 
  mot[which(mot == "NA")] <- NA
  
  #determine probabilities of person child with other people from population 
  res_ch <- child_pred(mother = mot, 
                             father = as.matrix(population[snps][1,]), snps = snps, trait = x, phenes = phen)
  #table with probabilities of possible children 
  comb_prob <- res_ch$comb_prob
  #table of probabilities of possible genotypes  
  comb_prob_2 <- res_ch$comb_prob_2
  
  #full tables
  for(i in 2:nrow(population)){
    res_ch <- child_pred(mother = mot, 
                father = as.matrix(population[snps][i,]), snps = snps, trait = x, phenes = phen)
    comb_prob <- rbind(comb_prob, res_ch$comb_prob)
    comb_prob_2 <- rbind(comb_prob_2, res_ch$comb_prob_2)
  }
  
  main_res <- stat_child(comb_prob = comb_prob, comb_prob_2 = comb_prob_2, phenes = phen, trait = x)
  main_res$probabilities <- main_res$probabilities/nrow(population)
  return(main_res)
})
names(result_ch) <- c("eyes","hair","skin")
result$child <- result_ch   

#make json-file
id=substr(basename(args$txt),1,6)
exportJson <- toJSON(result)
write(exportJson, file = paste0(args$out,"/",id,"_appearance.json"))
