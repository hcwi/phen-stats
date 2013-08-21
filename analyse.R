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
  
  if (file.exists("C:/strawberry/perl/bin/perl.exe")) {
    PERL <- "C:/strawberry/perl/bin/perl.exe"
  }
  
  studyAssayPairs <<- get.files()
  
  for (i in 1:dim(studyAssayPairs)) {
    
    sFile <- studyAssayPairs[i,1]
    aFile <- studyAssayPairs[i,2]
    
    experiment <<- load.files(sFile, aFile)   
    results <<- get.models(experiment)
    
    save.results(sFile, aFile, experiment, results)
    
  }
} 

# save results to files
save.results <- function (sFile, aFile, experiment, results) {
  
  print("[debug] save.results")
  
  sFile2 <- substr(sFile, start=0, stop=regexpr("[.]",sFile)-1)
  aFile2 <- substr(aFile, start=0, stop=regexpr("[.]",aFile)-1)
  
  rFile <- paste(sFile2, aFile2, "obj.R", sep="_")
  save(experiment, file=rFile)
  print(paste("R objects saved to file:", rFile))
  
  statFile <- paste(sFile2, aFile2, "stat.txt", sep="_")
  write.table(results, file=statFile)
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
  sa <- merge(study, assay, by=dupNames, all=T)
  print(paste("Merged by", dupNames))
  
  obs <- load.data(sa)
  d <<- obs
  dupNames <- subset(names(sa), match(names(sa), names(obs)) > 0)
  print(paste("Common columns: ", toString(dupNames)))
  sad <- merge(sa, obs, by=dupNames, all=T)
  print(paste("Merged by: ", toString(dupNames)))
  
  sad
}


get.models <- function(sad) {
  
  print("[debug] get.models")
  
  if(!require("lme4")) {
    install.packages("lme4")
  }
  library(lme4)
  
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
  
  
  results <- matrix(nrow=dim(effects$x)[2], ncol=0)
  
  models <- list()
  for (i in 1:length(traits)) {
    
    trait = traits[i] 
    print(paste("Models for a new trait:", trait))
    
    result <- matrix(ncol=2, nrow=0)
    colnames(result) <- c(paste(trait,"mean",sep=""), paste(trait,"s.e.", sep=""))
    
    form <- paste(trait,"~",fixed,"+(1|",random,")", sep="")
    print(paste("Formula:", form))
    
    model <- lmer(form, sadf)
    #model.coef <- coef(model)
    
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
      
      result <- rbind(result, cbind(means, sqrt(means.var)))
      
    }  
    
    #l <- list(trait=trait, fixed=fixed, random=random, model=model)
    #models <- c(models,list(l))
    results <- cbind(results, result)
  }
  
  results
  #models
}


prepare.effects <- function (sad) {
  
  print("[debug] prepare.effects")
  
  factorize <- function(x) {if (!is.numeric(x)) factor(x) else x}
  sadf <- lapply(sad, factorize)
  
  are.traits <- grep("Trait[.]Value", names(sad), value=T)
  
  are.levels <- grep("(Characteristics)|(Factor)", names(sad), value=T)
  have.var <- function(x) length(unique(x))>1
  are.var <- sapply(sad[are.levels], have.var)
  
  are.random <- grep("(Block)|(Field)|(Rank)|(Plot)", are.levels[are.var], value=T)
  are.fixed  <- grep("(Block)|(Field)|(Rank)|(Plot)", are.levels[are.var], value=T, invert=T)
  
  #analysis <- list(data=sadf, traits=are.traits, levels=are.levels, models)
  
  random <- are.random[1]
  if (length(are.random) > 1) {
    for (j in 2:length(are.random)) {
      random <- paste(random, are.random[j], sep="*")
    }
  }
  
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
  
  
  # generation of full model matrices for fixed effects
  print("Generation of full model matrices for fixed effects")
  
  # vector 1
  x.list <- list()
  x.full <- matrix(1, nrow=dim(sad)[1])
  colnames(x.full) <- "all"
  l.to <- 1
  x.list <- c(x.list, list(list(factor="const", from=1, to=l.to)))
  # single trait vectors
  for (j in 1:length(are.fixed)) {
    
    f <- paste(are.traits[1],"~",are.fixed[j],"+(1|",random,")-1", sep="")
    print(paste("  Tmp model formula:", f))
    mod <- lmer(f, sadf)
    m <- mod@X
    colnames(m) <- names(mod@fixef)
    x.full <- cbind(x.full, m)
    l.from <- l.to + 1
    l.to <- l.to + dim(m)[2]
    x.list <- c(x.list, list(list(factor=are.fixed[j], from=l.from, to=l.to)))   
  }  
  
  # all traits combination vector
  f <- paste(are.traits[1],"~",fixed.noconst,"+(1|",random,")-1", sep="")
  print(paste("  Tmp model formula:", f))
  mod <- lmer(f, sadf)
  m <- mod@X
  colnames(m) <- names(mod@fixef)
  x.full <- cbind(x.full, m)
  l.from <- l.to + 1
  l.to <- l.to + dim(m)[2]
  x.list <- c(x.list, list(list(factor=fixed.noconst, from=l.from, to=l.to)))
  
  x.full.unique <- unique(x.full)
  
  list(traits=are.traits, random=random, fixed=fixed, x=x.full.unique, x.list=x.list)
}



#Load data for study/assay pair
load.data <- function(sa) {
  
  print("[debug] load.data")
  
  dataName <- findDataFile(sa)
  d <- tryCatch({
    if (grep("xls", dataName)) {
      if (exists("PERL")) {
        print(paste("Loading", dataName,"with xls"))
        load.xls(paste(getwd(),dataName, sep="/"))
      }
      else {
        dataName2 <- gsub("([.]xls)|([.]xlsx)", ".txt", dataName)
        print(paste("Trying", dataName2,"with read.table"))
        read.table(dataName2, header=T, sep="\t")    
      }
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
  d <- read.xls(file, perl="C:/strawberry/perl/bin/perl.exe")
  d
}

#Find name of data file to use for study/assay pair
findDataFile <- function (sa) {
  
  print("[debug] findDataFile")
  
  are.files <- grep("(Raw)|(Derived)|(Processed).Data.File", names(sa), value=T)
  print(paste("Data files: ", toString(are.files)))
  
  have.all <- function(x) !any(is.na(x))
  are.full <- sapply(sa[are.files], have.all)
  
  if (length(are.full) < 1)
    stop("Among the columns referring to data there are no columns free of missing values")
  
  have.same <- function(x) length(unique(x))==1
  are.same <- sapply(sa[are.files][are.full], have.same)
  
  nSame = length(are.same)
  if ( nSame < 1)
    stop("Among the columns referring to data there are no full columns with all same values")
  if ( nSame > 1)
    warning("Among the columns referring to data there are more than one full columns with all same values. The last one will be used.")
  
  print(paste("Full equal data names in: ", toString(are.files[are.full][are.same])))
  dfName <- unique(sa[are.files][are.full][are.same][nSame])
  print(paste("Using data from file: ", dfName))
  dfName
}

#args <- commandArgs(TRUE)
#setwd(args[1])

options(stringsAsFactors=FALSE)
run()





calculate.means <- function(barley) {
  
  #UGLY: CALCULATIING MEANS BY HAND
  
  #select factors and phenotype by hand:
  #what <- "t.stemDiameter"
  #byWhat <- c("infraname", "f.block")
  
  #calculate means for all obs~factor combinations
  what <- names(barley[sapply(barley, is.numeric)])
  byWhat <- names(barley[sapply(barley, is.factor)])
  
  results <- list()
  for (j in 1:length(what)) {
    results[what[j]] <- mean(barley[,what[j]])
    for (i in 1:length(byWhat)) {
      mChar <- aggregate(barley[,what[j]]~barley[,byWhat[i]], FUN=mean)
      sdChar <- aggregate(barley[,what[j]]~barley[,byWhat[i]], FUN=sd)
      seChar <- aggregate(barley[,what[j]]~barley[,byWhat[i]], FUN=function(x) sd(x)/sqrt(length(x)))
      vChar <- cbind(mChar,sdChar[2],seChar[2])
      names(vChar) <- c(byWhat[i], paste(what[j], "mean", sep="."), paste(what[j],"sd", sep="."))
      results[[paste(what[j],byWhat[i],sep="~")]] <- vChar
    }
  }
  
  results
}

#P-VALS
#if(!require("languageR")) {install.packages("languageR")}
#library(languageR)
#m1.p <- pvals.fnc(m1)

#NULL MODEL COMPARISION
#m1.null <- lmer(t.length~1+(1|f.block), barley)
#anova(m1,m1.null)