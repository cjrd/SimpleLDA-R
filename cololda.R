# Script to run LDA
# first source auxfunctions.R
# then fill in inFile with the appropriate location of the input corpus data file
# the format of the input corpus data file is:
#  number_of_words_in_document1 word1:freq1 word2:freq2 etc
# where each line contains 1 document (see included corpus from Blei's LDA paper)
#
# after completion of LDA inference (takes several hours on computer with dual core processor and 4GB RAM)
# you will have an object 'ldamodel'
# str(ldamodel)
#
# Example to examine ldamodel:
# ldasumm = extract.topics(ldamodel, 20, 'data/vocab.txt') # third argument is the vocabulary map
# ldasumm$words
# ldasumm$probs
# ldasumm$logtopicprob


### SET THIS ###
inFile='data/ap.dat'


# # # initialization # # #
corpus = parse.data(inFile)

# important parameter 
k=10 # number of topics
# other parameters
estAlpha = TRUE
MAX.ES.IT = 20
ES.CONV.LIM = 1e-6
EM.CONV = 1e-3
MAX.EM = 50
alpha = 50/k;

# init the model randomly 
cwinit = matrix(runif(k*corpus$nterms),k,corpus$nterms) + 1/corpus$nterms
ldamodel = list(logProbW=matrix(rep(0,k*corpus$nterms),k,corpus$nterms), alpha = 1, ntopics=k, nterms=corpus$nterms)
sstats = list(ndocs=0,classword=cwinit,k,corpus$nterms,classtotal=rowSums(cwinit), alpha.ss = 0)
ldamodel = mstep.beta(ldamodel,sstats)
ldamodel$alpha = alpha 

like.hist = c() # keep track of the likelihood
likelihood.prev = 0
numit = 0
hasConverged = FALSE
nTopics = ldamodel$ntopics
nTerms = ldamodel$nterms
phi = matrix(rep(0,nTopics*nTerms),nTopics, nTerms)

# # # # Run variational expectation-maximization # # # #
while (!hasConverged){
  numit = numit + 1
  print(sprintf("----- EM Iteration %i ----- ", numit))
  
  # reset sufficient statistics and likelihood
  sstats$classword = matrix(rep(0,nTopics*nTerms), nrow=nTopics, nTerms) 
  sstats$classtotal = rep(0,nTopics)
  sstats$ndocs = 0
  sstats$alpha.ss = 0
  likelihood = 0
  
  # # # do E-step # # #
  for (d in 1:corpus$ndocs){
    if (d %% 200 == 0){
      print(sprintf("~~ completed e-step for %i docs ~~",d))
    } 
    
 # # do posterior inference # # 
    
    # initialize the document specific variables
    doc.oldlike = 0
    doc.length = corpus$docs[[d]]$dlength
    doc.totlen = corpus$docs[[d]]$total
    gammav = rep(ldamodel$alpha + doc.totlen/nTopics, nTopics)
    digamma.gam = rep(digamma(ldamodel$alpha + doc.totlen/nTopics), nTopics)
    phi = matrix(rep(1/nTopics, doc.length*nTopics), nrow=doc.length, ncol=nTopics)
    oldphi = phi[1,]
    
    # compute posterior dirichlet
    estep.converged = FALSE
    numits.es = 0;
    while (!estep.converged){
      numits.es = numits.es + 1
      # TODO: rewrite "R-style"
      for (n in 1:doc.length){
        phisum = 0
        for (k in 1:nTopics){
          oldphi[k] = phi[n,k]
          phi[n,k] = digamma.gam[k] + ldamodel$logProbW[k, corpus$docs[[d]]$words[[n]]+1]
          
          if (k > 1){
            phisum = log.sum(phisum,phi[n,k])
          }
          else{
            phisum = phi[n,k]
          }
        }
        
        for (k in 1:nTopics){
          phi[n,k] = exp(phi[n,k] - phisum)
          gammav[k] = gammav[k] + corpus$docs[[d]]$counts[[n]]*(phi[n,k] - oldphi[k])
          digamma.gam[k] = digamma(gammav[k])
          if (is.na(gammav[k])){
            print('error with gammav, contains na')
            browser()
          }
        }
      }
   
      # determine if the documents likelihood has converged
      doc.like = compute.likelihood(corpus$docs[[d]], ldamodel, phi, gammav)
      convfrac = (doc.oldlike - doc.like) / doc.oldlike
      doc.oldlike = doc.like


      if (convfrac < ES.CONV.LIM || numits.es > MAX.ES.IT){
        estep.converged = TRUE
        # print(sprintf("leaving E-step after %i iterations and convfrac: %1.3e, doc-likelihood: %1.3e", numits.es, convfrac, doc.like))
       # plot(doc.histlike)
      }
    } # end while e-step has not converged
    
  # # update the sufficient statistics for the M-step # #
    gamma.sum = sum(gammav)
    sstats$alpha.ss = sstats$alpha.ss + sum(sapply(gammav,digamma))
    sstats$alpha.ss = sstats$alpha.ss - nTopics*digamma(gamma.sum)
    
    for (n in 1:doc.length ){
      for (k in 1:nTopics){
        phink = corpus$docs[[d]]$counts[n]*phi[n,k]
        sstats$classword[k,corpus$docs[[d]]$words[n] + 1] = sstats$classword[k,corpus$docs[[d]]$words[n] + 1] + phink
        sstats$classtotal[k] = sstats$classtotal[k] + phink
      }
    }
    sstats$ndocs = sstats$ndocs + 1
    likelihood = likelihood + doc.like
  } # end for each document
  
# # # do M-step # # #
  
  print("[doing m-step]")
  
  # estimate beta
  ldamodel = mstep.beta(ldamodel,sstats)
  
  # estimate alpha
  if (estAlpha){
  D = sstats$ndocs
  alpha.iter = 0
  a.init = 100
  log.a = log(a.init)
  alpha.hasconv = FALSE
  while (!alpha.hasconv){
    alpha.iter = alpha.iter + 1
    a = exp(log.a)
    
    if (is.nan(a)){
      a.init = a.init*10
      print(sprintf("alpha became nan, initializing with alpha = %1.3e",a.init))
      a = a.init
      log.a = log(a)
    }
    
    f = D*(lgamma(nTopics*a) - nTopics*lgamma(a)) + (a-1)*sstats$alpha.ss
    df = D * (nTopics*digamma(nTopics*a) - nTopics*digamma(a)) + sstats$alpha.ss
    d2f = D * (nTopics*nTopics*trigamma(nTopics*a) - nTopics*trigamma(a))
    log.a = log.a - df/(d2f*a + df)
    print(sprintf("alpha optimization: %1.3e  %1.3e   %1.3e", exp(log.a), f, df))
    if (abs(df) < 1e-5 || alpha.iter > 100){
      alpha.hasconv = TRUE
    }
  }
  ldamodel$alpha = exp(log.a)
  }
  
  conv.em = (likelihood.prev - likelihood)/likelihood.prev
  likelihood.prev = likelihood
  like.hist[numit] = likelihood
  
  # make sure we're iterating enough for the likelihood to converge'
  if (conv.em < 0){ 
    MAX.ES.IT = MAX.ES.IT*2
  }
  if (((conv.em < EM.CONV && conv.em > 0)  || numit > MAX.EM) && numit > 2){
    print(sprintf("Converged with conv = %0.3f and %i iterations",conv.em,numit))
    hasConverged = TRUE
  }
  print(sprintf("likelihood: %1.4e, conv: %1.4e",likelihood, conv.em))
  plot(like.hist)
}
print("----- Finished -----")