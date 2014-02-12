# This file contains auxillary functions for the ColoLDA script

extract.topics <- function(ldamodel, nwords, dictfile){
# extract the top nwords words from each topic in beta using the dictionary provided in dictfile
#
# INPUT 
# ldamodel: the ldamodel list after using variational EM
# nwords: examine top nwords by likelihood
# dictfile: file containing the dictionary map of the word, each line contains one word, i.e.
#   line 1 contains 0->word, line 63 contains 64->word
#
# OUTPUT
# ldasumm: lda summary containing top nwords words for each topic, total log probabilities of each topic and
# probabilites of each of the top nwords for the topics
#
# author: colorado reed

indata = read.table(dictfile, stringsAsFactors = FALSE)
res = matrix(rep(0,ldamodel$ntopics*nwords), nrow=nwords, ncol=ldamodel$ntopics)
prbs = matrix(rep(0,ldamodel$ntopics*nwords), nrow=nwords, ncol=ldamodel$ntopics)
for (i in 1:ldamodel$ntopics){
  tmp = sort(ldamodel$logProbW[i,], decreasing=TRUE, index.return = TRUE)
  res[,i] = tmp$ix[1:nwords]
  prbs[,i] = tmp$x[1:nwords]
}
ldasumm = list()
ldasumm$words = matrix(indata[res], nrow=nwords, ncol=ldamodel$ntopics)
ldasumm$probs = exp(prbs)
ldasumm$logtopicprob = rowSums(ldamodel$logProbW)

return(ldasumm)
}

# # # log.sum # # #
log.sum <- function(log.a,log.b){
# sum of logarithms log(a+b)
  larger = max(c(log.a,log.b))
  smaller = min(c(log.a,log.b))
  res = larger + log(1 + exp(smaller - larger))
  return(res)
}

# # # compute.likelihood # # #
compute.likelihood <- function(doc, model, phi, gammav){
  # compute likelihood according to equation (15) in Blei's LDA paper
  likelihood = 0
  alpha = model$alpha
  nTopics = model$ntopics
  dig = sapply(gammav,digamma)
  gammav.sum = sum(gammav)
  digsum = digamma(gammav.sum)
  likelihood = lgamma(alpha*nTopics) - nTopics*lgamma(alpha) - lgamma(gammav.sum)
  # print(likelihood)
  for (k in 1:nTopics){
    addlike = (alpha - 1)*(dig[k] - digsum) + lgamma(gammav[k]) - (gammav[k] - 1)*(dig[k] - digsum)
    likelihood = likelihood 
    # print(sprintf("k_num %f",addlike))
    for (n in 1:doc$dlength){
      if (phi[n,k] > 0){
        addlike = doc$counts[n]*(phi[n,k]*((dig[k] - digsum) - log(phi[n,k]) + model$logProbW[k,doc$words[n]+1]))
        # print(sprintf("kn_num %f",addlike))
        likelihood = likelihood + addlike
      }
    }
  }
  # print(likelihood)
  return(likelihood)
}


mstep.beta <- function(ldamodel,sstats){
# estimate beta (logProbW) according to equation (7) of C.Reed's tutorial
  
  for (k in 1:ldamodel$ntopics){
    for (w in 1:ldamodel$nterms){
      if (sstats$classword[k,w] > 0 ){
        ldamodel$logProbW[k,w] = log(sstats$classword[k,w]) - log(sstats$classtotal[k])
      }
      else{
        ldamodel$logProbW[k,w] = -100
      }
    }
  }
  return(ldamodel)
}

parse.data <- function(inFile){
# parse document data for LDA
# author: colorado reed
  
  # init
  con  <- file(inFile, open = "r")
  num.docs <- 0
  docs <- list() #TODO better way to init?
  max.word <- -1
  
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    
    inVec <- (strsplit(oneLine, " "))
    
    if (length(inVec) == 0) next
    
    num.docs <- num.docs + 1
    len <- as.integer(inVec[[1]][1]) # number of terms
    words <- rep(0,len)
    wcounts <- rep(0,len)
    tot.words <- 0
    
    # acquire the term and term counts for each document
    for (i in 1:len){
      termvec <- strsplit(inVec[[1]][i+1],":");
      wcounts[i] <- as.integer(termvec[[1]][2])
      words[i] <- as.integer(termvec[[1]][1])
      
      if (words[i] > max.word){
        max.word <- words[i] # keep track of the largest term
      } 
      tot.words <- tot.words + wcounts[i]
    }
    docs[[num.docs]] <- list(words=words, counts=wcounts, dlength=len, total=tot.words)
    
  }
  
  close(con)
  corpus = list(docs=docs, nterms=max.word + 1, ndocs=num.docs) # +1 accounts for starting at zero
  return(corpus)
}