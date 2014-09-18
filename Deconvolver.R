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

options(digits= 22)
#scripts to compute pdfs from cdfs, do fourier magic and stuff...

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

getYfromfft <- function(fft, x) {
  x = x-1
  re = 0.0
  N = length(fft)
  for(t in 1:N) {
    re = re + fft[t] * exp((1i*2*pi*t*x)/N)
  }
  return( 1/N * re)
}


# markiere alle doppelten negativen Werte mit FALSE
negativeMirrors <- function(y){
  N = length(y)
  ret = rep(TRUE, N)
  
  thres = 1e-5
  
  for(i in 1:(N-1)){
    if(y[i] > 0 && y[i+1] < 0) {
      diff = y[i] + y[i+1]
      if(diff < thres){
        ret[i] = FALSE
        ret[i+1] = FALSE
      }
    }
  }
  return(ret)
}

domination <- function(y){
  N = length(y)
  range = 5
  
  d <- rep(0, N)
  
  for(i in 1:N){
    test = max(1, i-range):min(N, i+range)
    test = test[y[test]<y[i]]
    d[test] = d[test]+1
  }
  
  return(d<=range)
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


#returns the value of pdf at x plus the index where this value can be found.
getY <- function(pdf, x, min=1) {
  max = dim(pdf)[1]
  
  i = min
  while(i <= max) {
    if(pdf[i, 1] > x){
      return(c(pdf[i, 2], i))
    }
    i = i+1
  }
  return(c(pdf[max, 2], max))
  
}



makeSamples <- function(pdf, range=1e9) {
  cat("Sampling signal...\n")
  xs = seq(0, range, 1)
  
  re <- matrix(nrow=length(xs), ncol=2)
  
  idx = 1
  
  for(i in 1:length(xs)) {
    if(i %% 10000 == 0) {
      cat("\r", i/1e9 * 100, "%")
    }
    while(pdf[idx, 1] < xs[i]){
      idx=idx+1
    }
    re[i, 1] = xs[i]
    re[i, 2] = pdf[idx, 2]
    
  }
  cat("\n")
  
  return(re)
  
}

#pdf = f * c*g
#funktion gibt f aus (berechnet durch Faltungstheorem)
#pdf muss equally spaced sein. c gibt den faktor der beiden pdfs an.
unconvolve <- function(pdf, c=1) {
  
  #nullen anhaengen:
  zeros = matrix(nrow=(length(pdf)), ncol=1)
  for(i in 1:(length(pdf1))) {
    zeros[i] = 0
  }
  pdf = c(zeros, pdf, zeros)
  
  #berechne fourier transformation auf den samples:
  ff = fft(pdf)
  ff = sqrt(ff)/c     #wurzel ziehen
  ffti = abs(fft(ff, inverse=T) / length(ff))   #ruecktransformation
  
  #verschiebe ffti:
  ffti1 = ffti[1:(length(ffti)/2+1)]
  ffti2 = ffti[(length(ffti)/2+2):length(ffti)]
  ffti = c(ffti2, ffti1)
  
  return(ffti)
}


getSamples <- function(signal, numSamples) {
  stepsize = length(signal) / numSamples
  ret = matrix(nrow=numSamples, ncol=1)
  i = 1
  idx = 1
  while(i < length(signal)){
    #cat(i,"\n")
    cur = floor(i)
    nex = cur+1
    if(cur == nex){
      ret[idx] = signal[i]
      idx = idx+1
    } else {
      ab = 1-(i-cur)
      ret[idx] = ab * signal[cur] + (1-ab) * signal[nex]
      idx = idx+1
    }
    i = 1+(stepsize*(idx-1))
  }
  #cat(idx,"\n")
  ret[idx] = signal[length(signal)]
  return(ret[1:idx])
}



getAckPDFfromPay <- function(pay,c) {
  N = length(pay)
  
  ack = vector(length=N)
  
  for(n in N:1) {
    idx = ceiling(n*c)
    if(idx > N)
      next
    ack[idx] = ack[idx] + pay[n]
  }
  return(ack)
}



#z most NOT be a PDF but samples from that PDF!
gn <- function(t, z, beta) {
  N = length(z)
  s = 0+0i
  for(k in 1:N) {
    s = s + exp(1i*t*(z[k]/beta))
  }
  return(s/N)
}

#z is PDF
gn_pdf <- function(t, z, beta) {
  N = length(z)
  s = 0+0i
  #hist = c()
  hun = 0+0i
  
  #idx = 1:N
  idx = which(z>1e-20)
  cat(length(idx), "non-zero entries in gn_pdf...\n")
  i = 0
  for(k in idx) {
    s <- s + exp(1i*t*(k/beta)) * z[k]
    if(i %% 10 == 1) {
      suu = sum(abs(hun-s))
      if( suu < 1e-20) break;
      hun = s
    }
    #if(k %% 100 == 1) cat("gn_pdf", k, "/", N, s[100], z[k], suu, "\n")
    if(i %% 100 == 1) cat("gn_pdf", i, "/", length(idx), "\n")
    i=i+1
  }
  #plot(abs(hist))
  return(s)
}

f_dach_n <- function(t, gn_func, gamma, N) {
  cat("f_dach_n\n")
  oldP = 0+0i
  
  gn1 <- gn_func(gamma^(0) * t)
  gn2 <- gn_func(gamma^(1) * t)
  p <-  gn1 / gn2   #glied 0
  refID = 1
  for(k in 1:N) {
    oldP = p[refID]
    gn1 = gn_func(gamma^(2*k) * t)
    gn2 = gn_func(gamma^(2*k+1) * t)
    p <- p * ( gn1  / gn2 )
    if(abs(oldP - p[refID]) < 1e-40) break;
    cat("f_dach_n, ", k, "/", N, " - ", p[1], abs(oldP - p[1]), "\n")
  }
  cat("f_dach_n, ", k, "/", N, " - ", abs(oldP - p[1]), "\n")
  cat("\n")
  return(p)
}


characteristic_function_to_density <- function(
  phi, # characteristic function; should be vectorized
  n,   # Number of points, ideally a power of 2
  a, b # Evaluate the density on [a,b[
) {
  i <- 0:(n-1)            # Indices
  dx <- (b-a)/n           # Step size, for the density
  x <- a + i * dx         # Grid, for the density
  dt <- 2*pi / ( n * dx ) # Step size, frequency space
  c <- -n/2 * dt          # Evaluate the characteristic function on [c,d]
  d <-  n/2 * dt          # (center the interval on zero)
  t <- c + i * dt         # Grid, frequency space
  cat("evaluating density on [", a, ",", b, "[ or: ", c, ", ", d, "\n")
  phi_t <- phi(t)
  X <- exp( -(0+1i) * i * dt * a ) * phi_t
  Y <- fft(X)
  density <- dt / (2*pi) * exp( - (0+1i) * c * x ) * Y
  data.frame(
    i = i,
    t = t,
    characteristic_function = phi_t,
    x = x,
    density = Re(density)
  )
}



sampleCDF <- function(cdf, numSamples) {
  ret <- matrix(nrow=numSamples, ncol=2)
  
  stepSize = log10(length(cdf[,1]))/numSamples
  for(i in 1:numSamples){
    
    j = floor(10^(i*stepSize))
    ret[i,1] = cdf[j,1]
    ret[i,2] = cdf[j,2]
  }
  return(ret)
}



B_obs_intra2B_PL_intra <- function(){
  library(roxygen)
  cdf_B_obs_intra <- read.csv("B_obs_intra.csv", header = FALSE)
  pdf_B_obs_intra <- CDF2PDF(as.matrix(t(cdf_B_obs_intra)))    #input data
  
  
  #pdf_B_obs_intra <- getSimple_cdf_B_obs_intra()    #input data
  cdf_B_obs_intra <- PDF2CDF(pdf_B_obs_intra)
  beta <- 1
  gamma <- (66/(2.5*1448))
  
  
  maxFileSize = 1.25e9
  factor <- 100       #bin by 100 bytes
  vectorLength = maxFileSize / factor          #length of all bined vectors
  
  #make sure that vectorLength is good for convolving! for example 1e5 and 2.5e6 are bad but 3e6 is very good - has something todo with the fft implementation... takes years on bad values.
  
  pdfLong <- rep(0, vectorLength)
  
  
  for(i in 1:nrow(pdf_B_obs_intra)) {
    if(pdf_B_obs_intra[i, 2] > 0) {
      bytes <- pdf_B_obs_intra[i, 1]    #bytecount
      pr <- pdf_B_obs_intra[i, 2]       #probability of that bytecount
      
      if(bytes > vectorLength*factor)         #omit all bytecounts that are too large
        next
      
      idx <- ceiling(bytes/factor)            #bin index
      pdfLong[idx] <- pdfLong[idx] + pr
    }
  }
  
  pdfLong = c(pdfLong[1:vectorLength], rep(0, vectorLength))  #append zeros to pdf because this is the convolved vector
  
  pdfLong = pdfLong/sum(pdfLong)                              #normalize to get a pdf
  
  gn_z <- Curry(gn_pdf, z=pdfLong, beta=1)                    #setup g_n
  N <- length(pdfLong)                                        #N is now 2*vectorLength
  
  f <- Curry(f_dach_n, gn_func=gn_z, gamma=gamma, N=N)        #setup f_n
  d <- characteristic_function_to_density(f, N, 1, N)         #deconvolve the pdf
  
  x <- d$x
  y <- d$density
  
  
  manualCleaning = FALSE
  
  #########################################
  ######### manual cleaning phase #########
  #########################################
  if(manualCleaning) {
    
    t <- y[10:100]
    t[t<0.0004] = 0
    y[10:100] = t
    
    t <- y[200:1000]
    t[t<0.0014] = 0
    y[200:1000] = t
    
    t <- y[1000:1000000]
    t[t<0.0005] = 0
    y[1000:1000000] = t
    
    t <- y[1000000:10000000]
    t[t<2e-4] = 0
    y[1000000:10000000] = t
    
    #     idz = c(2,4,6,8,9,10,12,14,15,16,18,19)
    #     y[idz] = 0
    #     
    #     t <- y[21:75]
    #     t[t<0.01] = 0
    #     y[21:75] = t
    #     
    #     t <- y[20:500]
    #     t[t<0.005] = 0
    #     y[20:500] = t
    #     
    #     t <- y[76:3500]
    #     t[t<0.0005] = 0
    #     y[76:3500] = t
    #   
    #     t <- y[3501:10000]
    #     t[t<2e-4] = 0
    #     y[3501:10000] = t
    #     
    #     t <- y[10000:100000]
    #     t[t<3e-5] = 0
    #     y[10000:100000] = t
    #     
    #     t <- y[100000:300000]
    #     t[t<1e-5] = 0
    #     y[100000:300000] = t
    #     
    #     t <- y[30000:900000]
    #     t[t<1e-6] = 0
    #     y[30000:900000] = t
    #     
    #     t <- y[900000:1500000]
    #     t[t<2e-6] = 0
    #     y[900000:1500000] = t
    #     
    #     y[1500000:length(y)] = 0
    
    
    threshold = 1e-8
  } else {
    threshold = 1e-8                                          #everything below is treated as zeros
  }
  
  
  x = x[1:vectorLength] * factor                            #remove the zeros we appended before
  y = y[1:vectorLength]
  
  ys <- y[1:10]
  ys[ys<0] = 0
  y[11:length(y)] = y[11:length(y)] / (sum(y[11:length(y)]) / (1-sum(y[1:10])))
  
  
  
  
  dom = domination(y)
  neg = negativeMirrors(y)
  y[(!dom) | (!neg) | y<threshold] = 0
  y[y<threshold] = 0
  y[1:10] = ys
  y[11:length(y)] = y[11:length(y)] / (sum(y[11:length(y)]) / (1-sum(y[1:10])))
  
  #build pdf:
  pdfPAY = matrix(nrow=length(y), ncol=2)      #build the pdf-object - should be from 1 to vectorLength and in bins of size factor
  pdfPAY[,1] = x
  pdfPAY[,2] = y/sum(y)
  
  
  
  #build cdf:
  cdfPAY = PDF2CDF(pdfPAY)
  
  
  #create acks and reconstruct ALL traffic from previously extracted Payloads:
  #build cdf for ACK traffic:
  pdfACK = matrix(nrow=vectorLength, ncol=2)
  pdfACK[,1] = seq(factor, vectorLength*factor, factor)
  pdfACK[,2] = getAckPDFfromPay(pay=pdfPAY[,2], c=gamma)
  cdfACK = PDF2CDF(pdfACK)
  
  
  #build pdf for ALL traffic:
  cat("convolving two vectors of size", length(pdfPAY[,2]), "and", length(pdfACK[,2]), " - this may take some time...\n")
  c = convolve(pdfPAY[,2], rev(pdfACK[,2]), type="o")
  pdfALL = matrix(nrow=vectorLength, ncol=2)
  pdfALL[,1] = seq(factor, vectorLength*factor, factor)
  pdfALL[,2] = c[1:vectorLength]
  cdfALL = PDF2CDF(pdfALL)
  
  
  
  
  
  idx = c(seq(1, 1000, 1), seq(2000, 1e7, 1000))
  plot(cdf_B_obs_intra, log="x", xlim=c(1, 1.25e9), type="l", lwd=2, lty="solid", xlab="Flowsize", ylab="CDF")
  lines(cdfALL[idx,1], cdfALL[idx,2], col="red", lwd=2, lty="dotted")
  lines(cdfPAY[idx,1], cdfPAY[idx,2], col="blue", lwd=2, lty="dotdash")
  lines(cdfACK[idx,1], cdfACK[idx,2], col="green", lwd=2, lty="twodash")
  legend(1e6, 0.7, c("Original", "Payload", "Acks", "Convolved"),text.width=2, lty=c("solid", "dotted", "dotdash", "twodash"), col=c("black", "red", "blue", "green"), lwd=2)
  

  
}

B_obs_inter2B_PL_inter <- function(){
  library(roxygen)
  cdf_inRackFileSize <- read.csv("B_obs_inter.csv", header = FALSE)
  pdf_inRackFileSize <- CDF2PDF(as.matrix(t(cdf_inRackFileSize)))    #input data
  
  cdf_inRackFileSize <- PDF2CDF(pdf_inRackFileSize)
  beta <- 1
  gamma <- (66/(2.5*1448))
  
  
  maxFileSize = 1.25e9
  factor <- 100       #bin by 100 bytes
  vectorLength = maxFileSize / factor          #length of all bined vectors
  
  #make sure that vectorLength is good for convolving! for example 1e5 and 2.5e6 are bad but 3e6 is very good - has something todo with the fft implementation... takes years on bad values.
  
  pdfLong <- rep(0, vectorLength)
  
  for(i in 1:nrow(pdf_inRackFileSize)) {
    if(pdf_inRackFileSize[i, 2] > 0) {
      bytes <- pdf_inRackFileSize[i, 1]    #bytecount
      pr <- pdf_inRackFileSize[i, 2]       #probability of that bytecount
      
      if(bytes > vectorLength*factor)         #omit all bytecounts that are too large
        next
      
      idx <- ceiling(bytes/factor)            #bin index
      pdfLong[idx] <- pdfLong[idx] + pr
    }
  }
  
  pdfLong = c(pdfLong[1:vectorLength], rep(0, vectorLength))  #append zeros to pdf because this is the convolved vector
  
  pdfLong = pdfLong/sum(pdfLong)                              #normalize to get a pdf
  
  gn_z <- Curry(gn_pdf, z=pdfLong, beta=1)                    #setup g_n
  N <- length(pdfLong)                                        #N is now 2*vectorLength
  
  f <- Curry(f_dach_n, gn_func=gn_z, gamma=gamma, N=N)        #setup f_n
  d <- characteristic_function_to_density(f, N, 1, N)         #deconvolve the pdf
  
  x <- d$x
  y <- d$density
  
  
  
  manualCleaning = TRUE
  
  #########################################
  ######### manual cleaning phase #########
  #########################################
  if(manualCleaning) {
    
    
    
    t <- y[1e3:1e4]
    t[t<0.0008] = 0
    y[1e3:1e4] = t
    
    
    
    
    
    #      t <- y[1e5:1e6]
    #      t[t<0.0001] = 0
    #      y[1e5:1e6] = t
    
    #     t <- y[100000:200000]
    #     t[t<0.00022] = 0
    #     y[100000:200000] = t
    #     
    #     
    #     t <- y[300000:3000000]
    #     t[t<2e-4] = 0
    #     y[300000:3000000] = t
    #     
    
    #     idz = c(2,4,6,8,9,10,12,14,15,16,18,19)
    #     y[idz] = 0
    #     
    #     t <- y[21:75]
    #     t[t<0.01] = 0
    #     y[21:75] = t
    #     
    #     t <- y[20:500]
    #     t[t<0.005] = 0
    #     y[20:500] = t
    #     
    #     t <- y[76:3500]
    #     t[t<0.0005] = 0
    #     y[76:3500] = t
    #   
    #     t <- y[3501:10000]
    #     t[t<2e-4] = 0
    #     y[3501:10000] = t
    #     
    #     t <- y[10000:100000]
    #     t[t<3e-5] = 0
    #     y[10000:100000] = t
    #     
    #     t <- y[100000:300000]
    #     t[t<1e-5] = 0
    #     y[100000:300000] = t
    #     
    #     t <- y[30000:900000]
    #     t[t<1e-6] = 0
    #     y[30000:900000] = t
    #     
    #     t <- y[900000:1500000]
    #     t[t<2e-6] = 0
    #     y[900000:1500000] = t
    #     
    #     y[1500000:length(y)] = 0
    
    ####y[50:length(y)] = y[50:length(y)] / (sum(y[50:length(y)]) / (1-sum(y[1:50])))
    threshold = 1e-8
  } else {
    threshold = 1e-8                                          #everything below is treated as zeros
  }
  
  
  top = 12
  
  x = x[1:vectorLength] * factor                            #remove the zeros we appended before
  y = y[1:vectorLength]
  
  ys <- y[1:top] * 0.75
  ys[ys<0] = 0
  y[(top+1):length(y)] = y[(top+1):length(y)] / (sum(y[(top+1):length(y)]) / (1-sum(y[1:top])))
  
  
  dom = domination(y)
  neg = negativeMirrors(y)
  y[(!dom) | (!neg) | y<threshold] = 0
  y[y<threshold] = 0
  y[1:top] = ys
  y[(top+1):length(y)] = y[(top+1):length(y)] / (sum(y[(top+1):length(y)]) / (1-sum(y[1:top])))
  
  #build pdf:
  pdfPAY = matrix(nrow=length(y), ncol=2)      #build the pdf-object - should be from 1 to vectorLength and in bins of size factor
  pdfPAY[,1] = x
  pdfPAY[,2] = y/sum(y)
  
  #build cdf:
  cdfPAY = PDF2CDF(pdfPAY)
  
  
  #create acks and reconstruct ALL traffic from previously extracted Payloads:
  #build cdf for ACK traffic:
  pdfACK = matrix(nrow=vectorLength, ncol=2)
  pdfACK[,1] = seq(factor, vectorLength*factor, factor)
  pdfACK[,2] = getAckPDFfromPay(pay=pdfPAY[,2], c=gamma)
  cdfACK = PDF2CDF(pdfACK)
  
  
  #build pdf for ALL traffic:
  cat("convolving two vectors of size", length(pdfPAY[,2]), "and", length(pdfACK[,2]), " - this may take some time...\n")
  c = convolve(pdfPAY[,2], rev(pdfACK[,2]), type="o")
  pdfALL = matrix(nrow=vectorLength, ncol=2)
  pdfALL[,1] = seq(factor, vectorLength*factor, factor)
  pdfALL[,2] = c[1:vectorLength]
  cdfALL = PDF2CDF(pdfALL)
  
  
  idx = c(seq(1, 1000, 1), seq(2000, 1e7, 1000))
  plot(cdf_inRackFileSize, log="x", xlim=c(1, 1e12))
  lines(cdfALL[idx,1], cdfALL[idx,2], col="red")
  lines(cdfPAY[idx,1], cdfPAY[idx,2], col="blue")
  lines(cdfACK[idx,1], cdfACK[idx,2], col="green")
  
}