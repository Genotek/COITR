#===========================================================
# In:  SNPs of parents
# Out: Probabilities of trait variants and SNPs for children
#===========================================================
child_pred <- function(mother, father, snps, trait, phenes){
  #give names to parents SNPs
  names(mother) <- snps
  names(father) <- snps
  
  #the number of groups for prediction
  size <- length(names(phenes[[trait]][[1]][[1]]))
  
  #the list of probabilities of children's genotypes based on parents' genotypes
  genotypes <- vector(mode ="list", length =length(snps))
  names(genotypes) <- snps
  
  genotypes <-(lapply(snps, function(x){
    vsp <- vector(mode = "list", length = 2)
    names(vsp) <- c("parents", "child")
    
    #fill parents' genotypes
    vsp$parents <- vector(mode = "list", length = 2)
    names(vsp$parents) <- c("mother", "father")
    vsp$parents$father <- father[x]
    vsp$parents$mother <- mother[x]
    
    #make list for probabilities of children's genotypes
    vsp$child <- vector(mode = "list", length = 3)
    names(vsp$child) <- names(phenes[[trait]][x][[1]][1:3])
    
    #if one of the parents has "NA" genotype, child will have "NA" genotype too.
    if(sum(is.na(vsp$parents)) > 0){
      
      vsp$child[1:3] <- NA
      return(vsp)}
    
    else{
      
      #get parents alleles  
      mom <- str_split(vsp$parents$mother, pattern = "")[[1]]
      dad <- str_split(vsp$parents$father, pattern = "")[[1]]
      
      #all combinations of parents's alleles
      all <- as.vector(outer(mom, dad, paste0))
      
      #count frequancy of each combination
      all <- unlist(sapply(all, function(y){
        return(paste0(str_sort(str_split(y, pattern = "")[[1]]), collapse = ""))}))
      vsp$child <- as.list(sapply(names(vsp$child), function(z){
        a <- sum(all == z)/4
        return(a)}))
      
      return(vsp)}
  }))
  
  names(genotypes) <- snps
  
  for(i in 1:size){
    #filter snps with NA
    snps_clear <- snps[!is.na(mother)&!is.na(father)]
    
    #make list of genotype probabilities
    genes_frame <- vector(mode = "list", length = length(snps_clear))
    names(genes_frame) <- snps_clear
    
    #make list of conditional probabilities
    genes_frame_cond <- vector(mode = "list", length = length(snps_clear))
    names(genes_frame_cond) <- snps_clear
    
    for(x in snps_clear){
      v = rep(NA,3)
      for(j in 1:3){
        v[j] = phenes[[trait]][[x]][[j]][[i]]
      }
      genes_frame[[x]] <- unlist(genotypes[[x]]$child)
      genes_frame_cond[[x]] <- v
    }
    
    #In every iteration we make a vector, which consits of all multiplications of combinations
    #of snps, which number = number of iteration
    comb <- genes_frame[[1]][genes_frame[[1]] != 0]
    comb_2 <- genes_frame_cond[[1]][genes_frame[[1]] != 0]
    if(length(genes_frame) > 1){
      for(j in 2:length(genes_frame)){
        comb <- as.vector(outer(comb, genes_frame[[j]][genes_frame[[j]]!=0]))
        comb_2 <- as.vector(outer(comb_2, genes_frame_cond[[j]][genes_frame[[j]]!=0]))
      }
    }
    
    if(i == 1){
      #tables with probabilities
      comb_prob <- matrix(ncol = size, nrow = length(comb))
      comb_prob_2 <- matrix(ncol = size, nrow = length(comb))
    }
    
    comb_prob[,i] <- comb
    comb_prob_2[,i] <- comb_2
  }
  
  main_res <- vector(mode = "list", length = 2)
  names(main_res) <- c("comb_prob", "comb_prob_2")  
  
  main_res$comb_prob <- comb_prob
  main_res$comb_prob_2 <- comb_prob_2
  
  return(main_res)
}


#===========================================================
# In:  trait and SNP probabilities for children
# Out: mean, quantiles and SD of trait variant probabilities
#===========================================================
stat_child <- function(comb_prob_2, comb_prob, phenes, trait){
  #the number of groups for prediction
  size <- length(names(phenes[[trait]][[1]][[1]]))
  #the final list
  main_res <- vector(mode = "list", length = 3)
  names(main_res) <- c("probabilities", "quantiles", "sd")  
  
  #vector with probabilities of each possible phenotype
  res <- rep(0, size)
  for(i in 1:length(comb_prob[,1])){
    res <- res + (comb_prob_2[i,]*comb_prob[i,])/sum(comb_prob_2[i,])
    comb_prob_2[i,] <- comb_prob_2[i,]/sum(comb_prob_2[i,])
  }
  names(res) <- names(phenes[[trait]][[1]][[1]])
  main_res$probabilities <- res
  
  #the list with quantiles
  main_res$quantiles <- vector(mode = "list", length = size)
  names(main_res$quantiles) <- names(phenes[[trait]][[1]][[1]])  
  main_res$sd <- rep(NA, size)  
  names(main_res$sd) <- names(phenes[[trait]][[1]][[1]])
  for(i in 1:size){
    #20 - the number of quantiles, which are determined with the step 0.05
    main_res$quantiles[[i]] <- qdiscrete(c(0:20)*0.05, probs = comb_prob[,i][order(comb_prob_2[,i])], 
                                          values = sort(comb_prob_2[,i]))
    names(main_res$quantiles[[i]]) <- paste0(c(0:20)*5, "% quantile")
    main_res$sd[i] <- wt.sd(x = comb_prob_2[,i], wt = comb_prob[,i])
  }
  return(main_res)
}

#===========================================================
# In:  person's SNPs
# Out: person's SNPs, which we need
#===========================================================

read_snps <- function(chip, snps, chr_pos, gen){
  
  gen[snps] <- lapply(snps, function(x){
    
    #Find snps, which exist 
    if(sum(chip$id == x) > 0){ 
      cur_gen <- unique(chip$type[which(chip$id == x)])
      if(length(cur_gen) == 1){
        if(cur_gen == "--"){  #    
          return(NA)}         #
        else{                #
          return(as.character(cur_gen))}} 
      else{
        return(NA)}
    }
    else{
      return(NA)} 
  })
  
  #Try to get snps with the help of chrms and positions in it
  gen[snps][is.na(gen[snps])] <- lapply(snps[is.na(gen[snps])], function(x){
    if(sum(chip$chr == chr_pos[[x]]["chr"] & chip$pos == chr_pos[[x]]["pos"]) > 0){ 
      cur_gen <- unique(chip$type[which(chip$chr == 
                                          chr_pos[[x]]["chr"] & chip$pos == chr_pos[[x]]["pos"])])
      if(length(cur_gen) == 1){
        if(cur_gen == "--"){  #    
          return(NA)}         #
        else{                #
          return(as.character(cur_gen))}} 
      else{
        return(NA)}
    }
    else{
      return(NA)} 
  })
  return(gen)
}

#===========================================================
# In:  person's SNPs
# Out: probabilities of trait variants
#===========================================================
det_phen <- function(ph, chip){
  
  #This function will return this list 
  res <- vector(mode = "list", length = 3)
  names(res) <- c("Snps", "Probabilities", "Predictions") 
  
  #Get snps's names
  snps <- names(ph)[substring(names(ph), 1, 2) %in% "rs"] 
  genes <- names(ph)[substring(names(ph), 1, 3) %in% "gen"]
  
  #Get the human's snps
  gen <- rep(NA, length(c(snps, genes)))
  names(gen) <- c(snps, genes)
  gen <- read_snps(snps = snps, chr_pos = ph$chr_pos, chip = chip, gen = gen)
  
  #sorting genotype nucleotides in alphabetical order, e. g. TA will become AT after sorting
  gen <- unlist(gen)
  gen <- sapply(gen, function(x){
    paste0(str_sort(str_split(x, pattern = "")[[1]]), collapse = "")
  })
  
  #if phenotype needs in grouped genes
  if(length(genes) > 0){
    
    #Get the human's snps for haplotypes
    hapl <- rep(NA, 8)
    hapl <- read_snps(snps = c(ph$gen_MC1R_1$snps, ph$gen_MC1R_2$snps), 
                      chr_pos = ph$chr_pos, chip = chip, gen = hapl)
    
    # find the number of recessive alleles  
    m <- matrix(ncol = 8, nrow = 2)
    m[1,] <- unlist(hapl)
    m[2,] <- c("A", "A", "T", "T",  "C","T", "A", "A")
    vsp <- apply(m, 2, function(y){str_count(y[1],pattern = y[2]) })
    k <- c(NA, NA)
    o <- c(NA, NA)
    k[1] <- ifelse(sum(is.na(vsp[c(3:4)])) > 0, NA, 1)
    k[2] <- ifelse(sum(is.na(vsp[c(6:7)])) > 0, NA, 1)
    for(i in 1:2){
      if(sum(is.na(vsp[c(i*3, i*3 + 1)])) == 0){
        for(l in 1:4) {
          if(!is.na(vsp[l + 4*(i - 1)])){k[i] <- k[i] + vsp[l + 4*(i - 1)]}
        }  
        if(k[i] > 2){o[i] <- "RR"}
        if(k[i] == 2){o[i] <- "Rw"}
        if(k[i] == 1){o[i] <- "ww"}
      }
    }
    gen[genes] <- o
  }
  
  res$Snps <- gen
  
  # Check the presence of SNPs which are crucial for making the prediction
  # (if they are not, accuracy will degrade significantly)
  if(length(ph$main_rs) > 0){
    if(sum(is.na(gen[ph$main_rs])) > 0){
      res$Probabilities <- "We can't predict probabilities"
      res$Predictions <- "We can't predict colour"
      return(res)}
  }
  
  #Determine probabilities of each of the colour with NB model
  pr <- 1:length(ph[[1]][[1]])
  pr <- sapply(pr, function(y){ 
    prod(unlist(sapply(1:length(gen), function(x){
      return(ph[names(gen)[x]][[1]][gen[x]][[1]][[y]])
    })), na.rm = T)
  })
  pr <- pr/sum(pr)
  names(pr) <- names(ph[[1]][[1]])
  
  #For red colour
  if((length(genes) > 0) & (length(names(ph[[1]][[1]]))  >  2)){
    if(!is.na(gen[genes[1]])){
      if(gen[genes[1]] == "RR"){
        pr <- c(0, 0, 1, 0)
        names(pr) <- names(ph[[1]][[1]])
      }
    }
  }
  if(length(names(ph[[1]]))  >  2){
    res$Probabilities <- pr}
  
  #Use the threshold to determine, how many variants of phenotype are the most probable
  if(length(ph$threshold) > 0){
    if(max(pr) < ph$threshold){
      a <- pr
      a[which.max(pr)] <- 0
      res$Predictions <- c(names(pr)[which.max(pr)], names(pr)[which.max(a)])
      if (sum(c(names(pr)[which.max(pr)], names(pr)[which.max(a)]) %in% names(pr)[c(length(pr), 1)]) == 2){
        res$Predictions <- "We can't predict this phenotype"
      }
      return(res)
    }
    else{
      if(length(names(ph[[1]]))  >  2){
        res$Predictions <- names(pr)[which.max(pr)]
        return(res)}
    }
  }
  return(res)
}


drop.null <- function(x){
  if(is.null(x)) NA
  else x
}
