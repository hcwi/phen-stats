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
find.inv <- function(folder) {
  nums <- grep("^i_.*", list.files(path=folder))
  if (length(nums) == 0) {
    stop(paste("No investigation files were found in folder", folder))
  }
  if (length(nums) > 1) {
    stop(paste("More than one (", length(nums), ") investigation files were found in folder", folder))
  }
  inv <- list.files(path=folder)[nums]
  inv
}

#Run finding files
get.files <- function(dir=".") {
  inv <- read.inv(dir);
  files <- find.files(inv)
  files
}

#Run reading files # main
run <- function() {
  
  #setwd(paste(getwd(),"isatab", sep="/"))
  studyAssayPairs <<- get.files()
  
#   TODO processing and saving of each study+assay set
#   for (i in 1:dim(files)) {
#     load.files(files[i,1], files[i,2])
#     ..analyse..
#   }
  experiment <<- load.files(studyAssayPairs[1,1], studyAssayPairs[1,2])
  
  models <<- get.models(experiment)
}

#Load metadata for study/assay pair
load.files <- function(sName, aName) {
  
  print(paste(sName, aName))
  study <- read.table(sName, header=T)
  assay <- read.table(aName, header=T)
  
  dupNames <- subset(names(study), match(names(study), names(assay)) > 0)
  print(paste("Common columns: ", toString(dupNames)))
  dupNames.nontrivial <- grep("(REF)|(Accession.Number)", dupNames, invert=T)
  dupNames <- dupNames[dupNames.nontrivial]
  sa <- merge(study, assay, by=dupNames, all=T)
  print(paste("Merged by", dupNames))

  obs <- load.data(sa)
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
  
  factorize <- function(x) {if (!is.numeric(x)) factor(x) else x}
  sadf <- lapply(sad, factorize)
    
  are.traits <- grep("Trait[.]Value", names(sad), value=T)
  
  are.levels <- grep("(Characteristics)|(Factor)", names(sad), value=T)
  have.var <- function(x) length(unique(x))>1
  are.var <- sapply(sad[are.levels], have.var)
  
  are.random <- grep("(Block)|(Field)", are.levels[are.var], value=T)
  are.fixed  <- grep("(Block)|(Field)", are.levels[are.var], value=T, invert=T)
  
  models <- list()
  
  analysis <- list(data=sadf, traits=are.traits, levels=are.levels, models)
  
  for (i in 1:length(are.traits)) {
    
    trait = are.traits[i]
    
    fixed = are.fixed[1]
    for (j in 2:length(are.fixed)) {
      fixed <- paste(fixed, are.fixed[j], sep="*")
    }
    
    random <- are.random[1]
    if (length(are.random) > 1) {
      for (j in 2:length(are.random)) {
        random <- paste(random, are.random[j], sep="*")
      }
    }
    
    form <- paste(trait,"~",fixed,"+(1|",random,")", sep="")
    print(paste("Formula:", form))
  
    model <- lmer(form, sadf)
    model.coef <- coef(model)
    #print(model)
    
    
    l <- list(trait=trait, fixed=are.fixed, random=are.random, model=model)
    models <- c(models,list(l))
  }
  models

}



#Load data for study/assay pair
load.data <- function(sa) {
  
    dataName <- findDataFile(sa)
    d <- tryCatch({
      if (grep("xls", dataName)) {
        load.xls(paste(getwd(),dataName, sep="/"))
      }
      else {
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
  
  are.files <- grep("Data.File", names(sa), value=T)
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
  


#rm(list=ls())
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