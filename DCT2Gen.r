#
# DCT2Gen version 0.02
#
# For details, see the paper Philip Wette, Holger Karl, DCT2Gen: A Versatile TCP Traffic Generator for Data Centers, 2014.
#
# We highly recommend to carefully read the paper before using the traffic generator.
#
# Copyright (c) 2014 Philip Wette (wette@mail.upb.de)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# 
# If you need a different licensing, please contact Philip Wette (wette@mail.upb.de)
#
#


#invert a step function
reverseTable <- function(data) {
  tmp <- data[1, ]
  data[1, ] <- data[2, ];
  data[2, ] <- tmp;
  return(data)
}


# sample a value from a inverted cdf
getSample <- function(cdf_inv, scale = 1) {
  r <- runif(1, 0.0, 1.0)
  return(getVal(cdf_inv, r) * scale)
}


getVal <- function(ar, x) {
  #cat("getVal evaluating at ",x, "which is",floor(x*length(ar[1,]))+1, "and thus", ar[2, floor(x*length(ar[1,]))+1])
  return(ar[2, floor(x*length(ar[1,])) + 1  ]  )
}




#create traffic matrix using an ILP
createTM_ILP <- function(numRacks, hostsPerRack, trafficScale = 1, cdf_N_obs_intra_inv, cdf_N_obs_inter_inv, cdf_B_PL_intra_inv, cdf_B_PL_inter_inv) {
  servers <- numRacks * hostsPerRack;
  
  inRackPartners <- vector(length = servers)
  outOfRackPartners <- vector(length = servers)
  
  tm <- matrix(nrow = servers, ncol=servers);
  for(x in 1:servers) {
    #how many inRack und outRack communications does server i have?
    inRackPartners[x]    <- round(getSample(cdf_N_obs_intra_inv, hostsPerRack));
    outOfRackPartners[x] <- round(getSample(cdf_N_obs_inter_inv, 1500));    #TODO: since we do not know how things scale with increasing size of a data center, we assume it's ok to use the absolute numbers
    for(y in 1:servers) {
      tm[x,y] = 0
    }
  }
  
  cat("building a TM with", sum(inRackPartners), "in rack; and", sum(outOfRackPartners), "out rack TM entries...\n")
  
  #build a graph which has node degrees as stated by inRackPartners and outRackPartners:
  
  #in rack edges:
  
  #dump in rack degrees to file
  cat("", sep="", file="/tmp/ILPinputDegrees.csv", append=FALSE)
  for(x in 1:servers) {
    cat(inRackPartners[x], "\n", sep="", file="/tmp/ILPinputDegrees.csv", append=TRUE)
  }
  
  #run ILP:
  com <- paste("java -cp /opt/gurobi562/linux64/lib/gurobi.jar:./intraRackTMcreator.jar intraRackTMcreator", servers,  hostsPerRack, "/tmp/ILPinputDegrees.csv", "/tmp/tm.csv", sep=" ")
  system(com)
  Sys.sleep(5)
  
  #this is necessary beacuse of a bug that prevents R from reading a file where the first line does not contain any characters but \n
  system('echo XXX > /tmp/ilptmp; cat /tmp/tm.csv >> /tmp/ilptmp; rm /tmp/tm.csv; mv /tmp/ilptmp /tmp/tm.csv')
  
  #read result from file:
  ncol <- max(count.fields("/tmp/tm.csv", sep = ","));
  f <- read.table("/tmp/tm.csv", sep=",", header=T, blank.lines.skip=FALSE, fill=TRUE, col.names=paste0("V", seq_len(ncol)))
  
  ecIntra = 0
  
  #build traffic matrix:
  for(src in 1:servers) {
    for(dst in f[src,]) {
      if(is.na(dst) || dst < 0) {
        next
      }
      
      #we are building the intra rack matrix here - so obviously no inter rack edge must be produced. Let's check:
      if(!isInSameRack(src, dst+1, hostsPerRack)) {
        cat("ERROR!!! ILP computes bogus values!", src, dst+1, "\n")
      }
      
      tm[src, (dst+1)] = getSample(cdf_B_PL_intra_inv, trafficScale)
      tm[(dst+1), src] = getSample(cdf_B_PL_intra_inv, trafficScale)
      
      if(tm[src, (dst+1)] <= 0) {
        cat("Error: sampled zero tm entry size!\n")
      }
      if(tm[(dst+1),src] <= 0) {
        cat("Error: sampled zero tm entry size!\n")
      }
      
      ecIntra = ecIntra+2
    }
  }
  #that's it.
  cat("build", ecIntra, "intra rack edges.\n")
  
  
  numIn = 0
  numOut = 0
  #verify:
  for(src in 1:servers) {
    for(dst in 1:servers) {
      if(src != dst && tm[src,dst] > 0) {
        if(isInSameRack(src, dst, hostsPerRack)) {
          numIn = numIn+1
        } else {
          numOut = numOut+1
        }
      }
    }
  }
  cat("before intra shuffle: ", numIn, "in rack and", numOut, "out rack entries. total traffic:", sum(tm) ,"\n")
  
  #return(NULL)
  
  #shuffle edges:
  
  #iterate over all racks:
  for(x in 1:numRacks) {
    
    #used for shuffeling edges
    edges <- list()
    
    start = (x-1) * hostsPerRack + 1
    end = x * hostsPerRack
    
    cat("starting with rack", start, " to ", end, "\n")
    
    for(i in start:end){
      for(j in start:end){
        if(i == j) { next }
        if(i > j) { next }
        if(  tm[i,j] > 0 || tm[j,i] > 0) {
          edges[[length(edges)+1]] <- list("x" = i, "y" = j)
          #cat("+edge: ", i, j, "\n")
        }
      }
    }
    
    cat("shuffeling inRack edges...\n")
    #shuffle edges!
    for(s in 1:length(edges)) {
      e1 = sample.int(length(edges), 1)
      e2 = sample.int(length(edges), 1)
      tries = 0
      while(e1 == e2  || 
              edges[[e1]]$x == edges[[e2]]$y ||  #all four endpoints must be distinct
              edges[[e1]]$x == edges[[e2]]$x ||  #all four endpoints must be distinct
              edges[[e2]]$x == edges[[e1]]$y ||  #all four endpoints must be distinct
              edges[[e2]]$y == edges[[e1]]$y ||  #all four endpoints must be distinct
              tm[edges[[e1]]$x, edges[[e1]]$y] <= 1 || #and the selected edges must be present in the tm - TODO: this check is to hide a bug...
              tm[edges[[e2]]$x, edges[[e2]]$y] <= 1 || #and the selected edges must be present in the tm - TODO: this check is to hide a bug...
              tm[edges[[e1]]$x, edges[[e2]]$y] > 0 ||  #and the newly created edges must not already exist
              tm[edges[[e2]]$x, edges[[e1]]$y] > 0) {  #and the newly created edges must not already exist
        e2 = sample.int(length(edges), 1)
        
        tries = tries+1
        if(tries > length(edges)) {
          break
        }
      }
      if(tries > length(edges)) {
        next
      }
      
      # ei --w1--> ek
      # ei <-w2--- ek
      ei = edges[[e1]]$x
      ek = edges[[e1]]$y
      
      # ej --w3--> el
      # ej <-w4--- el
      ej = edges[[e2]]$x
      el = edges[[e2]]$y
      
      w1 = tm[ei,ek]
      w2 = tm[ek,ei]
      w3 = tm[ej,el]
      w4 = tm[el,ej]
      
      if(w1 <= 1 || w2 <= 1 || w3 <= 1 || w4 <= 1) {
        cat("Error! one of the w's is zero!", w1, w2, w3, w4, ei, ej, ek, el, "\n")
        return(NULL)
      }
      
      #delete old entries from tm:
      tm[ei,ek] = 0
      tm[ek,ei] = 0
      tm[ej,el] = 0
      tm[el,ej] = 0
      
      #create new edges:
      # ei --w1--> el
      # ei <-w4--- el
      # ej --w3--> ek
      # ej <-w2--- ek
      tm[ei,el] = w1
      tm[el,ei] = w4
      tm[ej,ek] = w3
      tm[ek,ej] = w2
      
      #cat("  adding (", ei, ",", el, ",", w1, ",", w4, ") and (", ej, ",", ek, ",", w3, ",", w2,")\n")
      
      #remove old edges from list:
      emax = e1
      emin = e2
      if(emin > emax){
        emax = e2
        emin = e1
      }
      #cat("-edge: ", edges[[max]]$x, edges[[max]]$y, "\n")
      #cat("-edge: ", edges[[min]]$x, edges[[min]]$y, "\n")
      
      edges[emax] <- NULL
      edges[emin] <- NULL
      
      
      #add new edges to list
      edges[[length(edges)+1]] <- list("x" = ei, "y" = el)
      edges[[length(edges)+1]] <- list("x" = ej, "y" = ek)
      
      #cat("+edge: ", edges[[length(edges)-1]]$x, edges[[length(edges)-1]]$y, "\n")
      #cat("+edge: ", edges[[length(edges)]]$x, edges[[length(edges)]]$y, "\n")
      
    }
    
  } #end for each rack.
  
  
  numIn = 0
  numOut = 0
  #verify:
  for(src in 1:servers) {
    for(dst in 1:servers) {
      if(src != dst && tm[src,dst] > 0) {
        if(isInSameRack(src, dst, hostsPerRack)) {
          numIn = numIn+1
        } else {
          numOut = numOut+1
        }
      }
    }
  }
  cat("after intra shuffle: ", numIn, "in rack and", numOut, "out rack entries. total traffic:", sum(tm) ,"\n")
  
  
  #out rack edges:
  edges <- list()
  
  #dump out rack degrees to file
  cat("", sep="", file="/tmp/ILPinputDegrees.csv", append=FALSE)
  for(x in 1:servers) {
    cat(outOfRackPartners[x], "\n", sep="", file="/tmp/ILPinputDegrees.csv", append=TRUE)
  }
  
  #run ILP:
  com <- paste("java -cp /opt/gurobi562/linux64/lib/gurobi.jar:./interRackTMcreator.jar interRackTMcreator", servers,  hostsPerRack, "/tmp/ILPinputDegrees.csv", "/tmp/tm.csv", sep=" ")
  system(com)
  Sys.sleep(5)
  
  #this is necessary beacuse of a bug that prevents R from reading a file where the first line does not contain any characters but \n
  system('echo XXX > /tmp/ilptmp; cat /tmp/tm.csv >> /tmp/ilptmp; rm /tmp/tm.csv; mv /tmp/ilptmp /tmp/tm.csv')
  
  #read result from file:
  ncol <- max(count.fields("/tmp/tm.csv", sep = ","));
  f <- read.table("/tmp/tm.csv", sep=",", header=T, blank.lines.skip=FALSE, fill=TRUE, col.names=paste0("V", seq_len(ncol)))
  
  ecInter = 0
  
  #build traffic matrix:
  for(src in 1:servers) {
    for(dst in f[src,]) {
      if(is.na(dst)) {
        next
      }
      
      #we are building the inter rack matrix here - so obviously no inter rack edge must be produced. Let's check:
      if(isInSameRack(src, dst+1, hostsPerRack)) {
        cat("ERROR!!! ILP computes bogus values!", src, dst+1, "\n")
      }
      
      tm[src, (dst+1)] = getSample(cdf_B_PL_inter_inv, trafficScale)
      tm[(dst+1), src] = getSample(cdf_B_PL_inter_inv, trafficScale)
      
      #if(tm[src, (dst+1)] <= 0) {
      #  cat("ERROR: inter rack traffic entry of size 0!", src, dst+1, "\n")
      #}
      #if(tm[(dst+1),src] <= 0) {
      #  cat("ERROR: inter rack traffic entry of size 0!", dst+1, src, "\n")
      #}
      
      edges[[length(edges)+1]] <- list("x" = src, "y" = (dst+1))
      
      ecInter = ecInter+2
    }
  }
  #that's it.
  
  cat("build", ecInter, "inter rack edges.\n")
  
  
  
  
  
  
  numIn = 0
  numOut = 0
  #verify:
  for(src in 1:servers) {
    for(dst in 1:servers) {
      if(src != dst && tm[src,dst] > 0) {
        if(isInSameRack(src, dst, hostsPerRack)) {
          numIn = numIn+1
        } else {
          numOut = numOut+1
        }
      }
    }
  }
  cat("before inter shuffle: ", numIn, "in rack and", numOut, "out rack entries. total traffic:", sum(tm) ,"\n")
  
  
  
  cat("shuffeling outRack edges...\n")
  #shuffle edges!
  for(s in 1:length(edges)) {
    
    cat(s , "/", length(edges), "\r")
    
    #oldTMsum = sum(tm)
    
    e1 = sample.int(length(edges), 1)
    e2 = sample.int(length(edges), 1)
    tries = 0
    while(e1 == e2  || 
            edges[[e1]]$x == edges[[e2]]$y ||  #all four endpoints must be distinct
            edges[[e1]]$x == edges[[e2]]$x ||  #all four endpoints must be distinct
            edges[[e2]]$x == edges[[e1]]$y ||  #all four endpoints must be distinct
            edges[[e2]]$y == edges[[e1]]$y ||  #all four endpoints must be distinct
            tm[edges[[e1]]$x, edges[[e2]]$y] > 0 ||  #and the newly created edges must not already exist
            tm[edges[[e2]]$x, edges[[e1]]$y] > 0 ||  #and the newly created edges must not already exist
            tm[edges[[e1]]$x, edges[[e1]]$y] <= 1 || #and the selected edges must be present in the tm - TODO: this check is to hide a bug...
            tm[edges[[e2]]$x, edges[[e2]]$y] <= 1 || #and the selected edges must be present in the tm - TODO: this check is to hide a bug...
            isInSameRack(edges[[e1]]$x, edges[[e2]]$y, hostsPerRack) || #and the newly created edges must not be intra rack edges
            isInSameRack(edges[[e1]]$y, edges[[e2]]$x, hostsPerRack) )  #and the newly created edges must not be intra rack edges
    {
      
      
      e2 = sample.int(length(edges), 1)
      
      tries = tries+1
      if(tries > 100) {
        break
      }
    }
    if(tries > 100) {
      #cat("not able to find a swapping edge - skipping.\n")
      next
    }
    
    # ei --w1--> ek
    # ei <-w2--- ek
    ei = edges[[e1]]$x
    ek = edges[[e1]]$y
    
    # ej --w3--> el
    # ej <-w4--- el
    ej = edges[[e2]]$x
    el = edges[[e2]]$y
    
    w1 = tm[ei,ek]
    w2 = tm[ek,ei]
    w3 = tm[ej,el]
    w4 = tm[el,ej]
    
    if(w1 <= 1 || w2 <= 1 || w3 <= 1 || w4 <= 1) {
      cat("Error! one of the w's is zero!", w1, w2, w3, w4, ei, ej, ek, el, "\n")
      return(NULL)
    }
    
    
    #delete old entries from tm:
    tm[ei,ek] = 0
    tm[ek,ei] = 0
    tm[ej,el] = 0
    tm[el,ej] = 0
    
    
    
    #create new edges:
    # ei --w1--> el
    # ei <-w4--- el
    # ej --w3--> ek
    # ej <-w2--- ek
    tm[ei,el] = w1
    tm[el,ei] = w4
    tm[ej,ek] = w3
    tm[ek,ej] = w2
    
    #remove old edges from list:
    emax = e1
    emin = e2
    if(emin > emax){
      emax = e2
      emin = e1
    }
    
    
    edges[emax] <- NULL
    edges[emin] <- NULL
    
    
    #add new edges to list
    edges[[length(edges)+1]] <- list("x" = ei, "y" = el)
    edges[[length(edges)+1]] <- list("x" = ej, "y" = ek)
    
    
    
  }
  cat("\nDone shuffeling inter rack edges.\n")
  
  
  numIn = 0
  numOut = 0
  #verify:
  for(src in 1:servers) {
    for(dst in 1:servers) {
      if(src != dst && tm[src,dst] > 0) {
        if(isInSameRack(src, dst, hostsPerRack)) {
          numIn = numIn+1
        } else {
          numOut = numOut+1
        }
      }
    }
  }
  cat("returning a TM with", numIn, "in rack and", numOut, "out rack entries. total traffic:", sum(tm) ,"\n")
  
  cat("wanted to build a TM with", sum(inRackPartners), "in rack and", sum(outOfRackPartners), "out rack TM entries...\n")
  
  
  
  return(tm);
}


isInSameRack <- function(x,y, hostsPerRack) {
  rx = floor((x-1) / hostsPerRack)
  ry = floor((y-1) / hostsPerRack)
  if(rx == ry)
    return(TRUE)
  
  return(FALSE)
}


createTMseries <- function(numMatrices = 6, numRacks = 75, hostsPerRack = 20, trafficScale = 1) {
  tms <- list()
  
  
  
  cat("reading N_obs_intra.csv\n")
  cdf_inRackComm <- read.csv("N_obs_intra.csv", header = FALSE)
  cdf_N_obs_intra_inv <- optimizeCDF(reverseTable(cdf_inRackComm))
  
  cat("reading N_obs_inter.csv\n")
  cdf_N_obs_inter <- read.csv("N_obs_inter.csv", header = FALSE)
  cdf_N_obs_inter_inv <- optimizeCDF(reverseTable(cdf_N_obs_inter))
  
  cat("reading B_PL_intra.csv\n")
  cdf_B_PL_intra <- read.csv("B_PL_intra.csv", header = FALSE)
  cdf_B_PL_intra_inv <- optimizeCDF(reverseTable(cdf_B_PL_intra))
  
  cat("reading B_PL_inter.csv\n")
  cdf_B_PL_inter <- read.csv("B_PL_inter.csv", header = FALSE)
  cdf_B_PL_inter_inv <- optimizeCDF(reverseTable(cdf_B_PL_inter))
  
  
  
  cat("done\n")
  
  
  
  for(i in 1:numMatrices) {
    tms[[i]] <- createTM_ILP(numRacks, hostsPerRack, trafficScale, cdf_N_obs_intra_inv, cdf_N_obs_inter_inv, cdf_B_PL_intra_inv, cdf_B_PL_inter_inv)
  }
  
  return(tms)
}

CDF2PDF <- function(cdf){
  
  c1 = c(cdf[,2], 1)
  c2 = c(0, cdf[,2])
  
  pdfy = c1 - c2
  pdfy = pdfy[1:(length(pdfy)-1)]
  
  pdf <- matrix(nrow=floor(dim(cdf)[1]), ncol=2)
  pdf[,2] = pdfy
  pdf[,1] = cdf[,1]
  
  return(pdf)
  
}

PDF2CDF <- function(pdf) {
  options(digits= 22)
  sum = 0.0
  
  cdf <- matrix(nrow=dim(pdf)[1], ncol=2)
  cdf[,1] = pdf[,1]
  cdf[1,2] = pdf[1,2]
  
  for(i in 2:dim(pdf)[1]) {    
    cdf[i,2] = cdf[i-1, 2] + pdf[i,2]
  }
  cdf[,2] = cdf[,2] / max(cdf[,2])
  
  return(cdf)
  
}




#create interarrival times + corresponding flow sizes
createInterTimesListWithTarget <- function(endtime, numServers, cdf_IAT_obs_inv, cdf_S_PL_inv, starttime, target){
  flows <- vector("list", numServers*500)	#TODO: just an educated guess of the size to prevent from subsequent reallocations
  i <- 1
  s <- 0
  
  stretchfactor <- (1500 / numServers)  # the cdf we use as input is ("normalized") for a cluster of 1500 servers.
  
  generatedBytes <- 0
  
  while(generatedBytes < target) {
    s <- s + ((getSample(cdf_IAT_obs_inv) / 1000.0) * stretchfactor)   #interarrival time is in milliseconds - devide by 1000 to convert to seconds.
    size <- getSample(cdf_S_PL_inv)
    flows[[i]] <- c(s, size)
    generatedBytes <- generatedBytes + size
    i <- i+1
  }
  i <- i-1
  lastTime <- flows[[i]][1]
  fak <- (endtime-starttime)/lastTime
  
  for(f in 1:i){
    flows[[f]][1] <- flows[[f]][1] * fak + starttime
  }
  
  #stretch flow start times such that it's in the interval starttime:endtime
  
  
  return(flows[1:i])
}



createFlowSet <- function(tm, numServers, cdf_IAT_obs_inv, cdf_S_PL_inv, round = 1) {
  
  cat("counting...")
  #count number of non zero entries & sum of traffic of the TM:
  nonZeros = 0;
  
  nonZeroE <- vector("list", numServers)
  
  trafficTM = 0
  for(x in 1:numServers) {
    nonZeroE[[x]] <- vector()
    for(y in 1:numServers) {
      if(tm[[x,y]] > 0) {
        nonZeros <- nonZeros+1
        trafficTM <- trafficTM + tm[[x,y]]
        nonZeroE[[x]] <- c(nonZeroE[[x]], y)
      }
    }
  }
  cat("\nTrafficTM (Payload only)", trafficTM, "\n")
  
  #create flows for the whole traffic Matrix:
  flows <- createInterTimesListWithTarget(10 * round, numServers, cdf_IAT_obs_inv, cdf_S_PL_inv, 10*(round-1), trafficTM);
  
  cat("done creating flows.\n")
  
  trafficFlows = 0
  for(i in 1:length(flows)) {
    trafficFlows <- trafficFlows + flows[[i]][2]
  }
  
  cat("Traffic Matrix:", trafficTM, "Traffic Flows:", trafficFlows, "Factor: ", trafficTM/trafficFlows, "\n")
  
  cat("We now have",nonZeros, "non zero entries in the traffic matrix and", length(flows) , "flows\n" )
  
  
  return(flows)
}


#tm is the traffic matrix, flows is NULL, numServers is the total number of servers, cdfs are cdfs, and round is the round (in 10 s)
mapFlowsToTM_DeficitRoundRobin <- function(tm, flows, numServers, cdf_IAT_obs_inv, cdf_S_obs_inv, round = 1) {
  returnFlows <- list()
  
  nonZeros = 0;
  
  nonZeroE <- vector("list", numServers)
  
  trafficTM = 0
  for(x in 1:numServers) {
    nonZeroE[[x]] <- vector()
    for(y in 1:numServers) {
      if(tm[[x,y]] > 0) {
        nonZeros <- nonZeros+1
        trafficTM <- trafficTM + tm[[x,y]]
        nonZeroE[[x]] <- c(nonZeroE[[x]], y)
      }
    }
  }
  cat("\nTrafficTM (Payload only)", trafficTM, "\n")
  
  #now we created flows for the traffic Matrix. We now have to assign these to the communication partners
  
  
  
  
  flowMatrix <- as.list(numeric(numServers^2))
  dim(flowMatrix) <- c(numServers,numServers)
  
  
  
  cat("starting RR phase...\n")
  #using round robin for flow assignment.
  quota <- matrix(ncol=numServers, nrow=numServers)
  got <- matrix(ncol=numServers, nrow=numServers)
  
  for(x in 1:numServers)
    for(y in 1:numServers) {
      flowMatrix[[x,y]] <- list()
      quota[[x,y]] <- 0
      got[[x,y]] <- 0
    }
  
  x = 1;
  y = 1;
  stat_round = 0;
  remainFak = 1;
  lastRoundF = 1
  toBigRefused <- 0
  
  yIdxInVec <- 1 # stores the index of the current y in nonZeroE[[x]]
  
  #test: sort flows:
  sizes = vector()
  for(f in 1:length(flows)) {
    sizes[f] = flows[[f]][2]
  }
  #sortedFlows = sort.list(sizes, decreasing=T)
  
  
  cat("writing debug output to file...\n")
  
  #write sizes to file
  write.csv(sizes, file="/tmp/flowSizes.csv")
  #write TM to file
  write.csv(tm, "/tmp/tm.csv", row.names=FALSE)
  
  cat("finding non-zero tm entries...\n")
  
  nonZeroPairsX <- vector()
  nonZeroPairsY <- vector()
  nonZeroPairsSize <- vector()
  nonzi = 1
  
  for(x in 1:numServers){
    for(y in 1:numServers){
      if(tm[x,y] > 0){
        nonZeroPairsX[nonzi] = x
        nonZeroPairsY[nonzi] = y
        nonZeroPairsSize[nonzi] = tm[x,y]
        nonzi = nonzi+1
      }
    }
  }
  
  
  #jeder tm entry bekommt initial einen moeglichst großen flow zugewiesen.
  if(FALSE) {
    cat("jedem eintrag einen flow geben. Wir haben initial", length(flows), "flows.\n")
    
    #sortiere tm eintraege absteigend
    sortedIds <- sort.list(nonZeroPairsSize, decreasing=T)
    cat("und", (nonzi-1), "nonzero tm eintraege.\n")
    
    #fuer jeden eintrag
    for(id in sortedIds){
      #suche naechst großen flow
      maxid = 1
      maxval = 0
      for(f in 1:length(flows)) {
        s = flows[[f]][2]
        if(s <= nonZeroPairsSize[id]){
          if(s > maxval){
            maxval = s
            maxid = f
          }
        }
      }
      if(maxval == 0)
        next
      
      x = nonZeroPairsX[id]
      y = nonZeroPairsY[id]
      #weise flow tm zu
      flowMatrix[[x,y]][[length(flowMatrix[[x,y]])+1]] <- flows[[maxid]]
      got[[x,y]] <- got[[x,y]] + flows[[maxid]][2]
      
      #loesche flow aus flows
      flows[[maxid]] <- NULL
    }
    
    cat("danach noch", length(flows), "uebrig\n")
  }  
  
  
  allFlows <- flows
  
  barriers <- c(1e11, 1e10, 1e9, 1e8, 1e7, 1e6, 1e5, 1e4, 1e3, 1e2, 1e1, 0)
  
  sizes = vector()
  for(f in 1:length(flows)) {
    sizes[f] = flows[[f]][2]
  }
  
  for(b in 1:(length(barriers)-1)){
    cat("starting with flows of size between", barriers[b], "and", barriers[b+1], "\n")
    indizes <- which(sizes<=barriers[b] & sizes>barriers[b+1])   #experimental
    
    if(length(indizes) < 1) {
      cat("no such flows.\n")
      next
    }
    
    flows <- allFlows[indizes]
    quota <- matrix(0, ncol=numServers, nrow=numServers)
    lastRoundF <- 1
    remainFak <- 1
    x = 1
    y = 1
    
    cat("which are", length(flows), "\n")
    
    for(f in 1:length(flows)) {
      #for(f in sortedFlows) {
      
      if(f %% 10000 == 1) {
        cat("\r", f, stat_round, toBigRefused)
      }
      found = FALSE
      while(!found) {
        remain <- tm[x,y] - got[[x,y]]
        quota[[x,y]] <- quota[[x,y]] +  max(0.01 * remain, 10)  #0.1,10 is a left shift of the distri. by approx 10% bytes
        if(quota[[x,y]] > flows[[f]][2] || remainFak > 1) {
          if( remain*remainFak > flows[[f]][2]) {
            flowMatrix[[x,y]][[length(flowMatrix[[x,y]])+1]] <- flows[[f]]
            quota[[x,y]] <- quota[[x,y]] - flows[[f]][2]
            got[[x,y]] <- got[[x,y]] + flows[[f]][2]
            
            if(remain < flows[[f]][2]){
              lost <- flows[[f]][2] - remain
              cat("loosing", lost, "bytes\n")
              
              # we will now remove the lost bytes from other tm entries of nearly the same size which has
              # still capacity left. We try to do that quite evenly accros all entries.
              # this is necessary because otherwise, if that happens with large entries, the remaining entries eat up all samll flows
              # and thus not enough small flows end up to be distributed to the small tm entries.
              
              # get capacity class of the tm entry
              capacityclass <- 1e11
              while(tm[x,y] < capacityclass) {
                capacityclass <- capacityclass*0.1
              }
              upCap <- capacityclass*10
              loCap <- upCap*0.1
              
              #find out all tm entries between upCap and loCap:
              idxx <- which(nonZeroPairsSize>=loCap & nonZeroPairsSize<upCap)
              totN <- 0
              for(d in idxx){
                loX <- nonZeroPairsX[d]
                loY <- nonZeroPairsY[d]
                if(tm[loX,loY] > got[[loX,loY]]){
                  
                  totN <- totN+1
                }
              }
              #we now know that we have totN entries left in this capacityclass. we will now take some of these bytes away.
              share <- lost / totN
              
              #cat("   we are in capacity class", capacityclass, "with", totN, "non full tm entries (out of",length(idxx),") -> distributing a share of", share, "\n")
              
              verschnitt <- 0
              totN <- 0
              tot <- 0
              for(d in idxx){
                loX <- nonZeroPairsX[d]
                loY <- nonZeroPairsY[d]
                if(tm[loX,loY] > got[[loX,loY]]){
                  thisshare <- share
                  remain <- tm[loX,loY] - got[[loX,loY]]
                  if(remain < share){
                    thisshare <- remain
                    verschnitt <- verschnitt + share - remain
                  } else {
                    totN <- totN+1
                  }
                  
                  #cat("    lowering", loX, ",", loY, "from", tm[loX, loY])
                  tm[loX, loY] <- tm[loX, loY] - thisshare
                  #cat("to", tm[loX, loY], "\n")
                  tot <- tot+remain
                }
              }
              
              #verschnitt fair aufteilen - je nachdem wer wie viel noch frei hat.
              if(verschnitt > 0) {
                for(d in idxx){
                  loX <- nonZeroPairsX[d]
                  loY <- nonZeroPairsY[d]
                  if(tm[loX,loY] > got[[loX,loY]]){
                    share <- (tm[loX,loY] - got[[loX,loY]])/tot
                    thisshare <- share*verschnitt
                    #cat("    lowering", loX, ",", loY, "from", tm[loX, loY])
                    tm[loX, loY] <- tm[loX, loY] - thisshare
                    #cat("to", tm[loX, loY], "\n")
                  }
                }
              }
              
              
              
            }
            
            remainFak <- 1
            found = TRUE
          } else {
            toBigRefused <- toBigRefused +1
          }
        }
        
        
        #next - only iterate over nonZero entries.
        #when a flow did iterate over all pairs without getting assigned, we calculate remainFak which is the 
        #factor of the largest entry in the residual matrix and the flowSize.
        
        #y <- min(vec[vec>y])                                          # high run time! 23%
        yIdxInVec <- yIdxInVec + 1
        y <- nonZeroE[[x]][yIdxInVec]
        if(is.na(y)) {
          while(is.na(y)) {
            x <- x+1
            if(x > numServers) {
              x <- 1
              stat_round <- stat_round +1
              #cat("we made another round. current flow", f, "lastRoundF", lastRoundF, "\n")
              
              if(f == lastRoundF) {
                remainM <- tm - got
                remainMax <- max(remainM)*remainFak
                
                if(remainMax <= 0) {
                  return(flowMatrix)
                }
                
                #if(remainMax < flows[[f]][2]) {
                remainFak <- flows[[f]][2] / remainMax +1
                lostBytes = max(flows[[f]][2] - remainMax, 0)
                cat(f, stat_round, toBigRefused, " Setting remainFak to", remainFak, "\n");
                #}
              }
              
              lastRoundF <- f
            }
            y <- nonZeroE[[x]][1]
            yIdxInVec <- 1
          }
          
        }
        
        
        
      }
      
    }
    cat("\n")
    
  }
  
  return(flowMatrix)
  
}



optimizeCDF <- function(cdf) {
  cat("optimizing table...\n")
  cdf[1, dim(cdf)[2]] = 1                 #set last probability value to 1
  idx <- 2;
  opt <- matrix(nrow=2, ncol=dim(cdf)[2])
  opt[1,1] <- cdf[1,1]
  opt[2,1] <- cdf[2,1]
  for(i in 2:dim(cdf)[2]) {                   #unify prop. values
    if(cdf[1, i-1] != cdf[1, i]) {      
      opt[1, idx] <- cdf[1, i]
      opt[2, idx] <- cdf[2, i]
      idx <- idx+1
    } else {
      opt[2, idx] <- cdf[2, i]
    }
  }
  
  opt2 <- matrix(nrow=2, ncol=idx)
  
  for(i in 1:idx-1) {
    opt2[1, i] <- opt[1, i]
    opt2[2, i] <- opt[2, i]
  }
  
  
  
  l <- 100000
  len <- dim(cdf)[2]
  opt3 <- matrix(nrow=2, ncol=l)
  
  
  p <- 1
  for(i in 1:l) {
    x <- i/l
    
    
    while(opt2[1, p] < x && p < dim(cdf)[2]) {
      p <- p+1
    }
    
    opt3[1,i] <- x;
    opt3[2,i] <- opt2[2,p];
  }
  
  
  #linearize the outcome:
  last = 1
  for(new in 2:l) {
    if(opt3[2,last] != opt3[2,new] && opt3[2,last] > 0) {
      
      c = new - last
      if(c > 1) {
        diff = opt3[2,new] - opt3[2,last]
        for(i in last:new-1) {
          opt3[2,i] = opt3[2,i] + diff * ((i - last)/c)
        }
      }
      
      last = new
    }
  }
  
  
  return(opt3)
}


writeFlowMatrixToFile <- function(matrix, fileName) {
  cat("writing flows to file...\n");
  for(x in 1:dim(matrix)[1]) {
    cat("\r", x, "... ")
    for(y in 1:dim(matrix)[2]) {
      if(length(matrix[[x,y]]) > 0) {
        for(idx in 1:length(matrix[[x,y]])) {
          f <- matrix[[x,y]][[idx]]
          cat(x,", ",y,", ",f[1],", ",f[2],"\n", sep="", file=fileName, append=TRUE)
        }
      }
    }
  }
  cat("\ndone\n")
}


makeTraffic <- function(numMatrices = 1, numRacks = 75, hostsPerRack = 20, trafficScale = 1, fileName = "/tmp/flows.csv") {
  tms <- createTMseries(numMatrices, numRacks, hostsPerRack, trafficScale);
  
  cat("reading S_obs.csv\n")
  cdf_S_obs <- read.csv("S_obs.csv", header = FALSE)
  cdf_S_obs_inv <- optimizeCDF(reverseTable(cdf_S_obs))
  
  cdf_S_PL <- cdf_S_obs2cdf_S_PL(cdf_S_obs_inv)
  cdf_S_PL_inv <- optimizeCDF(reverseTable(cdf_S_PL))
  
  cat("reading IAT_obs.csv\n")
  cdf_IAT_obs <- read.csv("IAT_obs.csv", header = FALSE)
  cdf_IAT_obs_inv <- optimizeCDF(reverseTable(cdf_IAT_obs))
  
  #create empty flow file:
  cat("", sep="", file=fileName, append=FALSE)
  
  for(round in 1:numMatrices) {
    flows <- createFlowSet(tms[[round]], numRacks * hostsPerRack, cdf_IAT_obs_inv, cdf_S_PL_inv, round)
    flowMatrix <- mapFlowsToTM_DeficitRoundRobin(tms[[round]], flows, numRacks * hostsPerRack, cdf_IAT_obs_inv, cdf_S_PL_inv, round)
    writeFlowMatrixToFile(flowMatrix, fileName);
  }
  
  
}

cdf_S_obs2cdf_S_PL <- function(cdf_inv){
  pdf <- CDF2PDF(t(reverseTable(cdf_inv)))
  x <- pdf[,1]
  y <- pdf[,2]
  
  N = length(y)
  
  n = N
  
  higherCandidates <- vector()
  nextJ <- -1
  
  oldAckSize <- -1
  ackSize <- -1
  
  cdfNotValid = F
  
  for(i in N:1) {
    paySize <- x[i]
    
    oldAckSize <- ackSize
    ackSize <- (ceiling(paySize / (2*1448))+4)*66 * 1.1
    
    # falls nicht mehr genug wahrscheinlichkeit bei ackSize uebrig ist, ziehen wir etwas von kleineren flows ab
    # diese duerfen aber nicht mehr als 30% kleiner sein.
    minAck <- 0
    remainProp <- y[i]
    
    if(remainProp <= 0)
      next
    
    
    #cat("decrement", paySize, ackSize, remainProp, "\n")
    
    if(ackSize > paySize){
      cat("Warning: some payload flows create larger ack flows than the payload flows itself. Not able to handle this situation.\n")
      cat("However, this operation still worked but the resulting flow size CDF may be biased!\n")
      break
    }
    
    
    while(remainProp > 0) {
      p = 0
      for(j in n:1){
        if(j <= 0) break;
        if(y[j] <= 0) next
        if(x[j] > ackSize) next
        
        
        if(nextJ < 0 | oldAckSize != ackSize)
          nextJ <- j+1
        
        #found possible candidate: find next higher one:
        ol <- nextJ
        for(l in nextJ:N) {
          if(y[l] > 0) {
            nextJ <- l
            break
          }
        }
        #if(i %% 1000 == 1) cat("took", (ol - nextJ), "inner, and ", (n-j), "outer loops. NextJ = ", nextJ, "\n")
        
        
        n = j
        
        #if(i %% 1000 == 1) cat("choose between", j, "&", nextJ, "meaning", x[j], "and", x[nextJ], "for ackSize", ackSize)
        p <- j
        if(abs(x[j] - ackSize) < abs(x[nextJ] - ackSize)) {
          #j
          p <- j
        } else {
          #nextJ
          p <- nextJ
        }
        
        #if(i %% 1000 == 1) cat(" - choosing", p, "\n")
        
        
        break;
        
        
      }
      
      if(p <= 0) {
        #we don't have a flow that is lower than ackSize. Die.
        break;
      }
      
      # found a possible candidate with non-zero prob
      if(y[p] >= remainProp) {
        y[p] <- y[p] - remainProp
        remainProp <- 0
      } else {
        remainProp <- remainProp - y[p]
        y[p] <- 0
      }
      
    }
    
    
    if(remainProp > 0) {
      #cat("missing a flow of size", ackSize, "to kill for the ack flow from", paySize, "\n")
      cdfNotValid = T
    }
    
  }
  
  pay <- matrix(nrow=N, ncol=2)
  pay[,1] <- x
  pay[,2] <- y/sum(y)
  
  if(cdfNotValid) {
    cat("The S_obs did not contain enough small flows. It was not possible to find S_PL that results in the given S_obs.\n")
    cat("This is just a warning. We carry on with a skewd payload size distribution.\n")
  }
  
  return(t(PDF2CDF(pay)))
  
}

cat("To compute a traffic schedule, use the function makeTraffic(..)")
