
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

args <- commandArgs(TRUE)
analyse(args[1])
