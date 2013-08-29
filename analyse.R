#Find pairs of study/assay file names in investigation file
find.files <- function(inv) {
  
  print("[debug] find.files")
  
  files <- data.frame(studyName=character(0), assayName=character(0), stringsAsFactors=FALSE)
  
  inv.studies <- as.matrix(subset(inv, V1=="Study File Name")) 
  inv.assays <- as.matrix(subset(inv, V1=="Study Assay File Name"))
  
  for (i in 1:dim(inv.studies)[1]) {
    sName <- inv.studies[i,2]
    for (j in 2:dim(inv.assays)[2]) {
      aName <- inv.assays[i,j]
      if (is.null(aName) || aName =="") {
        break
      }
      else {
        pair <- c(sName, aName)
        files <- rbind(files, pair)
      }
    }
  }
  files  
}

#Read standard investigation file in given folder
read.inv <- function(folder=".", file="i_Investigation.txt") {
  
  print("[debug] read.inv")
  
  inv = tryCatch({
    fname = paste(folder, file, sep="/")
    read.delim(fname, header=F)
  },    
                 error = function(e)       
                   print(e),
                 warning = function(w) {
                   f <- find.inv(folder)
                   inv=read.inv(folder, f)
                 },
                 finally = function() {
                   on.exit(close(fname))
                 }
  )
}

#Find non-standard investigation file in given folder
find.inv <- function(dir) {
  
  print("[debug] find.inv")
  
  nums <- grep("^i_.*", list.files(path=dir))
  if (length(nums) == 0) {
    stop(paste("No investigation files were found in folder", dir))
  }
  if (length(nums) > 1) {
    stop(paste("More than one (", length(nums), ") investigation files were found in folder", dir))
  }
  inv <- list.files(path=dir)[nums]
  inv
}

#Run finding files - first investigation, then study/assay pairs
get.files <- function(dir=".") {
  
  print("[debug] get.files")
  
  inv <- read.inv(dir);
  files <- find.files(inv)
  files
}

#Run processing: find, read, model, save
run <- function() {
  
  print("[debug] run")
  
  if (file.exists("C:/strawberry/perl/bin/perl.exe")) 
    PERL <<- "C:/strawberry/perl/bin/perl.exe"
  
  studyAssayPairs <<- get.files()
  
  for (i in 1:dim(studyAssayPairs)[1]) {
    
    sFile <- studyAssayPairs[i,1]
    aFile <- studyAssayPairs[i,2]  
    sa <- load.files(sFile, aFile)
    
    dFile <- find.dFile(sa)
    studyAssayPairs[i,3] <<- dFile
    dat <- load.data(dFile)
    d <<- dat[[1]]
    d.names <<- dat[[2]]
    
    experiment <<- get.experiment(sa, d)
    
    result <<- get.models(experiment, d.names)
    
    save.results(sFile, aFile, experiment, result)
    
  }
} 

# save results to files
save.results <- function (sFile, aFile, experiment, res) {
  
  print("[debug] save.results")
  
  means <- res[[1]]
  models <- res[[2]]
  
  sFile2 <- substr(sFile, start=0, stop=regexpr("[.]",sFile)-1)
  aFile2 <- substr(aFile, start=0, stop=regexpr("[.]",aFile)-1)
  
  rFile <- paste(sFile2, aFile2, "obj.R", sep="_")
  
  
  save(experiment, file=rFile)
  save(models, file=rFile)
  print(paste("R objects saved to file:", rFile))
  
  a <- read.table(aFile, header=T, check.names=F, sep="\t")
  names <- names(a)
  print(names)
  
  statFile <- paste(sFile2, aFile2, "stat.txt", sep="_")
  write.table(means, file=statFile, sep="\t", na="", row.names=F)
  print(paste("Sufficient statistics saved to file: ", statFile))
  
  # update isa-tab file to include sufficient data file
  update.file(aFile, statFile)   
  print(paste("Assay file", aFile, "updated to include sufficient statistics column"))
} 

update.file <- function(aFile, statFile) {
  
  print("[debug] update.file")
  
  print("Updating..")
  a <- read.table(aFile, header=T, check.names=F, sep="\t")
  a <- cbind(a, "Sufficient Data File"=statFile)
  write.table(a, na="", row.names=F, sep="\t", file=paste(aFile, 2, ".txt", sep=""))
  #write.table(a, na="", row.names=F, sep="\t", aFile)
  
}


#Load metadata for study/assay pair
load.files <- function(sName, aName) {
  
  print("[debug] load.files")
  
  print(paste(sName, aName))
  study <- read.table(sName, header=T, sep='\t')
  print("[debug]     read.table study")
  assay <- read.table(aName, header=T, sep='\t', fill=T)
  print("[debug]     read.table assay")
  s <<- study
  a <<- assay
  
  dupNames <- subset(names(study), match(names(study), names(assay)) > 0)
  print(paste("Common columns: ", toString(dupNames)))
  dupNames.nontrivial <- grep("(REF)|(Accession.Number)", dupNames, invert=T)
  dupNames <- dupNames[dupNames.nontrivial]
  warning("Parameter 'all' in merging changed to FALSE -- rows not matching study/assay will be removed. Rethink!")
  sa <- merge(study, assay, by=dupNames, all=F)
  print(paste("Merged by", dupNames))
  
  sa  
}


get.experiment <- function(sa, d) {
  
  print("[debug] get.experiment")
  
  dupNames <- subset(names(sa), match(names(sa), names(d)) > 0)
  print(paste("Common columns: ", toString(dupNames)))
  sad <- merge(sa, d, by=dupNames, all=T)
  print(paste("Merged by: ", toString(dupNames)))
  
  sad
}


get.models <- function(sad, d.names) {
  
  print("[debug] get.models")
  
  if(!require("lme4")) {
    install.packages("lme4")
  }
  library(lme4)
  if(!require("reshape")) {
    install.packages("reshape")
  }
  library(reshape)
  
  
  effects <- prepare.effects(sad)
  traits <- effects$traits
  random <- effects$random
  fixed <- effects$fixed
  
  
  
  #   factorize <- function(x) {if (!is.numeric(x)) factor(x) else x}
  #   sadf <- lapply(sad, factorize)
  
  factorizenames <- function(x) {
    if (is.numeric(sad[,x]) && 
          length(grep("Trait[.]Value", names(sad)[x])) > 0) 
      sad[,x] 
    else factor(sad[,x])
  }
  sadf <- lapply(seq_along(sad), factorizenames)
  names(sadf) <- names(sad)
  
  
  #factors <- strsplit(effects$x.list[[4]]$factor, "[*]")
  c_fixed <- unlist(strsplit(fixed, "[*]"))
  c_random <- unlist(strsplit(random, "[*]"))
  c_est <- paste("", unlist(strsplit(traits, "[*]")), sep="")
  c_se <- paste("S.e.", unlist(strsplit(traits, "[*]")), sep="")
  c_traits <- c(rbind(c_est,c_se))
  c_type <- "Parameter"
  cols <- c(c_type, c_fixed, c_random, c_traits)
  
  results <- matrix(nrow=dim(effects$x)[2], ncol=length(cols))
  colnames(results) <- cols
  
  results[,c_type] <- "Mean"
  
  
  for (j in 1:length(effects$x.list)) {
    
    factors <- strsplit(effects$x.list[[j]]$factor, "[*]")[[1]]
    from <- effects$x.list[[j]]$from
    to <- effects$x.list[[j]]$to
    
    if (factors[1] != "const") {
      for (i in from:to) {
        col <- strsplit(colnames(effects$x)[i], "_")[[1]]
        for (k in 1:length(col)) {
          results[i,factors[k]] <- col[k]
        }
      }
    }
  }
  
  models <- list()
  for (i in 1:length(traits)) {
    
    trait = traits[i] 
    print(paste("Models for a new trait:", trait))
    
    form <- paste(trait,"~",fixed,"+(1|",random,")", sep="")
    print(paste("Formula:", form))
    
    model <- lmer(form, sadf)
    
    x <- unique(model@X)
    est <- x %*% model@fixef
    
    est.cov <- x %*% vcov(model) %*% t(x)
    
    for (j in 1:length(effects$x.list)) {
      factor <- effects$x.list[[j]]$factor
      from <- effects$x.list[[j]]$from
      to <- effects$x.list[[j]]$to
      
      xf <- effects$x[,from:to]
      
      m <- solve(t(xf) %*% xf) %*% t(xf)
      means <- m %*% est
      rownames(means) <- colnames(effects$x)[from:to]
      
      means.var <- diag(m %*% est.cov %*% t(m))
      
      results[from:to, trait] <- means
      results[from:to, paste("S.e.",trait, sep="")] <- sqrt(means.var)
      
      
    }  
    
    l <- list(trait=trait, fixed=fixed, random=random, model=model)
    models <- c(models,list(l))
  }
  
  list(results, models)
}


prepare.effects <- function (sad) {
  
  print("[debug] prepare.effects")
  
  factorize <- function(x) {if (!is.numeric(x)) factor(x) else x}
  sadf <- lapply(sad, factorize)
  
  are.traits <- grep("Trait[.]Value", names(sad), value=T)
  
  warning("Removing traits with no variation")
  have.var <- function(x) length(unique(x))>1
  are.var <- sapply(sad[are.traits], have.var)
  are.traits <- are.traits[are.var]
  
  are.levels <- grep("((Characteristics)|(Factor))", names(sad), value=T)
  
  # Only for Keygene data testing
  warning("Filtering factors to exclude *id* names -- only for Keygene data. Remove for other analyses!")
  are.levels <- grep("[Ii]d", are.levels, value=T, invert=T)
  
  are.var <- sapply(sad[are.levels], have.var)
  
  are.random <- grep("(Block)|(Field)|(Rank)|(Plot)|(Replic)|(Column)|(Row)", are.levels[are.var], value=T)
  are.fixed  <- grep("(Block)|(Field)|(Rank)|(Plot)|(Replic)|(Column)|(Row)", are.levels[are.var], value=T, invert=T)
  
  if (length(are.random) > 0) {
    random <- are.random[1]
    if (length(are.random) > 1) {
      for (j in 2:length(are.random)) {
        random <- paste(random, are.random[j], sep="*")
      }
    }
  }
  else {
    random <- grep("Source.Name", names(sad), value=T)
  }
  print(paste("[debug]      random: ", random, sep=""))
  
  fixed = are.fixed[1]
  fixed.noconst = are.fixed[1]
  if (length(are.fixed) > 1) {
    for (j in 2:length(are.fixed)) {
      #if (j>2) {
      #  print(" -------- More than 2 fixed effects! System is not preapred to deal with multiple fixed factors which usually lead to  not positive definite matrices. Only 2 first factors will be used -------- ")
      #} else 
{
  fixed <- paste(fixed, are.fixed[j], sep="*")  
  fixed.noconst <- paste(fixed.noconst, are.fixed[j], sep=":")  
}
    }
  }
  print(paste("[debug]      fixed: ", fixed, sep=""))
  
  # pivot
  
  names <- vector()
  for (i in 1:length(are.fixed)) {
    
    com <- combn(are.fixed, i)
    for (j in 1:dim(com)[2]) {
      name <- com[1,j]
      if (dim(com)[1] > 1) {
        for (k in 2:dim(com)[1]) {
          name <- paste(name, com[k,j], sep="*")
        }
      }
      names <- rbind(names, name)
    }
    
  }
  
  
  x.list <- list()
  x.full <- matrix(1, nrow=dim(sad)[1])
  colnames(x.full) <- "all"
  l.to <- 1
  x.list <- c(x.list, list(list(factor="const", from=1, to=l.to)))
  
  for (i in 1:length(names)) {
    
    formula <- paste("Sample.Name", sep="~", names[i])
    print(paste("[debug]           names: ", names[i], sep=""))
    print(paste("[debug]           formu: ", formula, sep=""))
    
    x <- cast(sad, formula, length)
    x <- x[-1]
    
    x.full <- cbind(x.full, x)
    l.from <- l.to + 1
    l.to <- l.to + dim(x)[2]
    x.list <- c(x.list, list(list(factor=names[i], from=l.from, to=l.to)))   
  }
  
  x.full.unique <- as.matrix(unique(x.full))  
  
  list(traits=are.traits, random=random, fixed=fixed, x=x.full.unique, x.list=x.list)
}



#Load data for study/assay pair
load.data <- function(dFile) {
  
  print("[debug] load.data")
  print(paste("Loading", dFile))
  
  f <- paste(getwd(),dFile, sep="/")
  
  d <- tryCatch({
    if (grep("xls", dFile)) {
      if (exists("PERL"))
        load.xls(f)
      else
        load.txt(f)
    }
    else {
      print(paste("Loading", dataName,"with read.table"))
      read.table(dataName, header=T, sep="\t")
    }
  },    
                error = function(e)       
                  print(e),
                warning = function(w) {
                  print(w)
                },
                finally = function() {
                  on.exit(close(dataName))
                }
  )
  d
}

#Load data from xls file
load.xls <- function(file) {
  
  print("[debug] load.xls")
  
  if(!require("gdata")) {install.packages("gdata")}
  library(gdata)
  d <- read.xls(file, perl=PERL)
  d2 <- read.xls(file, perl=PERL, check.names=F)
  d.names <- as.vector(names(d2))
  names(d.names) <- names(d)
  list(d, d.names)
}

#Load data from txt file
load.txt <- function(file) {
  
  print("[debug] load.txt")
  
  dataName2 <- gsub("([.]xls)|([.]xlsx)", ".txt", dataName)
  print(paste("Trying", dataName2,"with read.table"))
  d <- read.table(dataName2, header=T, sep="\t")    
  d2 <- read.table(dataName2, header=T, sep="\t", check.names=F)    
  d.names <- as.vector(names(d2))
  names(d.names) <- names(d)
  list(d, d.names)
}


#Find name of data file to use for study/assay pair
find.dFile <- function (sa) {
  
  print("[debug] findDataFile")
  
  are.files <- grep("(Raw)|(Derived)|(Processed).Data.File", names(sa), value=T)
  print(paste("Data files: ", toString(are.files)))
  
  have.all <- function(x) !any(is.na(x))
  are.full <- sapply(sa[are.files], have.all)
  
  if (length(are.full) < 1)
    stop("Among the columns referring to data there are no columns free of missing values")
  
  have.same <- function(x) length(unique(x))==1
  are.same <- sapply(sa[are.files][are.full], have.same)
  
  nSame = length(sa[are.files][are.full][are.same])
  if ( nSame < 1)
    stop("Among the columns referring to data there are no full columns with all same values")
  if ( nSame > 1)
    warning("Among the columns referring to data there are more than one full columns with all same values. The last one will be used.")
  
  print(paste("Full equal data names in: ", toString(are.files[are.full][are.same])))
  dfName <- unique(sa[are.files][are.full][are.same][nSame])
  print(paste("Using data from file: ", dfName))
  dfName
}



# Things to do before running in Java

# remove global variables (<<-)
# change assay file instead of adding "*2.txt"

# uncomment:
# args <- commandArgs(TRUE)
# setwd(args[1])

options(stringsAsFactors=FALSE)
run()

#setwd("C:/Users/hcwi/Desktop/phen-stats/isatab")
#setwd("C:/Users/hcwi/Desktop/phen/src/test/resources/DataWUR")
#setwd("C:/Users/hcwi/Desktop/phen/src/test/resources/Phenotyping2")
#setwd("C:/Users/hcwi/Desktop/phen/src/test/resources/IPGPASData")
#setwd()

# calculate means for all obs~factor combinations
#   what <- names(barley[sapply(barley, is.numeric)])
#   byWhat <- names(barley[sapply(barley, is.factor)])
#   mChar <- aggregate(barley[,what[j]]~barley[,byWhat[i]], FUN=mean)
# 

#P-VALS
#if(!require("languageR")) {install.packages("languageR")}
#library(languageR)
#m1.p <- pvals.fnc(m1)

#NULL MODEL COMPARISION
#m1.null <- lmer(t.length~1+(1|f.block), barley)
#anova(m1,m1.null)