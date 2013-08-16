#Find pairs of study/assay file names in investigation file
find.files <- function(inv) {
  
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
  inv <- read.inv(dir);
  files <- find.files(inv)
  files
}

#Run processing: find, read, model, save
run <- function() {
  
  if (file.exists("C:/strawberry/perl/bin/perl.exe")) {
    PERL <- "C:/strawberry/perl/bin/perl.exe"
  }
  
  studyAssayPairs <<- get.files()
  
  for (i in 1:dim(studyAssayPairs)) {
    
    sFile <- studyAssayPairs[i,1]
    aFile <- studyAssayPairs[i,2]
    
    experiment <<- load.files(sFile, aFile)   
    models <<- get.models(experiment)
    
    save.results(sFile, aFile, experiment, models)
    
  }
} 

# save results to files
save.results <- function (sFile, aFile, experiment, models) {
  
  sFile2 <- substr(sFile, start=0, stop=regexpr("[.]",sFile)-1)
  aFile2 <- substr(aFile, start=0, stop=regexpr("[.]",aFile)-1)
  
#   outDir <- "output"
#   if (!file.exists(outDir)) {
#     print("Creating outDir")
#     dir.create(outDir)
#   }
#   rFile <- paste(outDir, sep="/", paste(sFile2, aFile2, "obj.R", sep="_"))
  rFile <- paste(sFile2, aFile2, "obj.R", sep="_")
  save(experiment, models, file=rFile)
  print(paste("R objects saved to file:", rFile))
  
#   statFile <- paste(outDir, sep="/", paste(sFile2, aFile2, "stat.txt", sep="_"))    
  statFile <- paste(sFile2, aFile2, "stat.txt", sep="_")
  for (i in 1:length(models)) {
    if(i==1)
      write.table(models[[i]]$model@fixef, file=statFile)
    else 
      write.table(models[[i]]$model@fixef, file=statFile, append=T)
  }
  print(paste("Sufficient statistics saved to file: ", statFile))
  
  # update isa-tab file to include sufficient data file
  update.file(aFile, statFile)   
  print(paste("Assay file", aFile, "udpdated to include sufficient statistics column"))
} 

update.file <- function(aFile, statFile) {
  print("Updating..")
  a <- read.table(aFile, header=T, check.names=F, sep="\t")
  a <- cbind(a, "Sufficient Data File"=statFile)
  write.table(a, na="", row.names=F, sep="\t", file=paste(aFile, 2, ".txt", sep=""))
  #write.table(a, na="", row.names=F, sep="\t", aFile)
  
}


#Load metadata for study/assay pair
load.files <- function(sName, aName) {
  
  print(paste(sName, aName))
  study <- read.table(sName, header=T)
  assay <- read.table(aName, header=T, fill=T, sep="\t")
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
  
  if(!require("lme4")) {
    install.packages("lme4")
  }
  library(lme4)
  
  
  effects <- prepareEffects(sad)
  random <- effects$random
  fixed <- effects$fixed
  x.full.u <- effects$x
  x.list <- effects$x.list
  
  factorize <- function(x) {if (!is.numeric(x)) factor(x) else x}
  sadf <- lapply(sad, factorize)
  
  models <- list()
  for (i in 1:length(are.traits)) {
    
    trait = are.traits[i] 
    
    form <- paste(trait,"~",fixed,"+(1|",random,")", sep="")
    print(paste("Formula:", form))
    
    model <- lmer(form, sadf)
    model.coef <- coef(model)

    
    
    l <- list(trait=trait, fixed=fixed, random=random, model=model)
    models <- c(models,list(l))
  }
  
  models
}


prepareEffects <- function (sadf) {
  
  factorize <- function(x) {if (!is.numeric(x)) factor(x) else x}
  sadf <- lapply(sad, factorize)
  
  are.traits <- grep("Trait[.]Value", names(sad), value=T)
  
  are.levels <- grep("(Characteristics)|(Factor)", names(sad), value=T)
  have.var <- function(x) length(unique(x))>1
  are.var <- sapply(sad[are.levels], have.var)
  
  are.random <- grep("(Block)|(Field)", are.levels[are.var], value=T)
  are.fixed  <- grep("(Block)|(Field)", are.levels[are.var], value=T, invert=T)
  
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
      fixed <- paste(fixed, are.fixed[j], sep="*")  
      fixed.noconst <- paste(fixed.noconst, are.fixed[j], sep=":")  
    }
  }
  
  
  # generation of full model matrices for fixed effects
  print("Generation of full model matrices for fixed effects")
  
  # vector 1
  x.list <- list()
  x.full <- matrix(1, nrow=dim(sad)[1])
  l.to <- 1
  x.list <- c(x.list, list(list(from=1, to=l.to)))
  
  # single trait vectors
  for (j in 1:length(are.fixed)) {
    
    f <- paste(are.traits[1],"~",are.fixed[j],"+(1|",random,")-1", sep="")
    print(paste("  Tmp model formula:", f))
    m <- lmer(f, sadf)@X
    x.full <- cbind(x.full, m)
    l.from <- l.to + 1
    l.to <- l.to + dim(m)[2]
    x.list <- c(x.list, list(list(trait=are.fixed[j], from=l.from, to=l.to)))   
  }  
  
  # all traits combination vector
  f <- paste(are.traits[1],"~",fixed.noconst,"+(1|",random,")-1", sep="")
  print(paste("  Tmp model formula:", f))
  m <- lmer(f, sadf)@X
  x.full <- cbind(x.full, m)
  l.from <- l.to + 1
  l.to <- l.to + dim(m)[2]
  x.list <- c(x.list, list(list(trait=fixed.noconst, from=l.from, to=l.to)))
  
  
  x.full.unique <- unique(x.full)
  
  list(random=random, fixed=fixed, x=x.full.unique, x.list=x.list)
}



#Load data for study/assay pair
load.data <- function(sa) {
  
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
  
  if(!require("gdata")) {install.packages("gdata")}
  library(gdata)
  d <- read.xls(file, perl="C:/strawberry/perl/bin/perl.exe")
  d
}

#Find name of data file to use for study/assay pair
findDataFile <- function (sa) {
  
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
























analyse <- function (path) {
  
  setwd(path)
  print(path)
  #TODO return errors when libraries are missing (they do not install on their own unless cran mirror is chosen)
  barley <- load.data()
  model <- model.data(barley)
  means <- calculate.means(barley)
  save(barley, model, means, file="output/savedObjects.R")
  fname = "output/stats.txt"
  write.table(model@fixef, file=fname)
  fname
}

load.data <- function () {
  
  #LOADING DATA
  s <- read.table("s_Study1.txt", header=T)
  a <- read.table("a_study1_phenotyping.txt", header=T)
  
  if(!require("gdata")) {install.packages("gdata")}
  library(gdata)
  
  # d <- read.xls("a_study1_processed_data.xlsx", perl="C:/strawberry/perl/bin/perl.exe")
  # alternative to gdata (which requires perl) is transformation to txt file:
  d <- read.table("a_study1_processed_data.txt", header=T, sep="\t")
  
  sa <- merge(s,a, by="Source.Name", all=T)
  sad <- merge(sa, d, by="Sample.Name", all=T)
  
  sad$Term.Accession.Number <- factor(sad$Term.Accession.Number)
  sad$Term.Accession.Number.1 <- factor(sad$Term.Accession.Number.1)
  sad$Term.Accession.Number.2 <- factor(sad$Term.Accession.Number.2)
  sad$Term.Accession.Number.3 <- factor(sad$Term.Accession.Number.3)
  sad$Term.Accession.Number.y <- factor(sad$Term.Accession.Number.y)
  sad$Term.Source.REF.3 <- factor(sad$Term.Source.REF.3)
  sad$Term.Source.REF.y <- factor(sad$Term.Source.REF.y)
  sad$Factor.Value.Block. <- factor(sad$Factor.Value.Block.)
  sad$Raw.Data.File <- factor(sad$Raw.Data.File)
  
  sad$Trait.value.Stem.diameter. <- levels(sad$Trait.value.Stem.diameter.)[as.integer(sad$Trait.value.Stem.diameter.)]
  
  barley.full <- sad
  barley <- barley.full[c(1,2,6,14,17,24,25,28)]
  names(barley) <- c("sample", "source", "infraname", "f.treatment", "f.block", "t.length", "t.colour", "t.stemDiameter")
  
  barley
}

model.data <- function(barley) {
  
  #CONSTRUCTING MODEL
  
  if(!require("lme4")) {
    install.packages("lme4")
  }
  library(lme4)
  model <- lmer(t.length~infraname*f.treatment+(1|f.block), barley)
  #model.coef <- coef(model)
  
  model
}

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

#CLEANUP
#rm("a", "s", "d", "sa", "sad", "i", "j", "what", "byWhat", "seChar", "sdChar", "vChar", "mChar")
#ls()

#RESULTS
#results
#model
#model.coef


#ADDITIONAL 

#P-VALS
#if(!require("languageR")) {install.packages("languageR")}
#library(languageR)
#m1.p <- pvals.fnc(m1)

#NULL MODEL COMPARISION
#m1.null <- lmer(t.length~1+(1|f.block), barley)
#anova(m1,m1.null)

#write(c(1:10),file='/output/stats.txt')
#f <- paste(args[1], '/output/stats.txt', sep='')
#write(c(1:10),file=f)

#args <- commandArgs(TRUE)
#analyse(args[1])