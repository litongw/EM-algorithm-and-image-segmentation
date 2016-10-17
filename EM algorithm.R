# set up working directory
setwd("~/R files")
# read histograms data
H<-matrix(readBin("histograms.bin", "double", 640000), 40000, 16)
dim(H) #results in a 40000*16 matrix, each row is a histogram with 16 bins
# Step 1.Implement the EM algorithm
MultinomialEM<-function(H,K,tau){
  if(0 %in% H){
    H=H+0.01
  }
  set.seed(3)
  n<-dim(H)[1]
  d<-dim(H)[2]
  # randomly choose K of histograms
  theta<-H[sample(1:n,K),]
  
  init<-1
  counter<-0
  a<-matrix(0,n,K)
  delta<-Inf
  
  while(init>=tau){
    # E-step
    a_old<-a
    phi<-exp(H %*% t(log(theta)))
    c<-rep(1/K,K)
    a<-t(t(phi)*c)
    a<-a/rowSums(a)
    a[is.na(a)]<-0
    # M-step
    b<-t(a)%*%H
    theta<-b/rowSums(b)
    if(init==1){
      init=0
      next
    }
    # Compute a measure of the change of assignments during the current iteration
    
    counter<-counter+1
    delta<-norm(a-a_old,type="1")
  }
  m<-apply(a,1,FUN = which.max)
  return(m)
}

# Step 2.Run the algorithm and prepare visulization fuction
visu<-function(H,K,tau){
  result<-MultinomialEM(H,K,tau)
  l<-sqrt(length(result))
  out<-matrix(result,l,l)[,200:1]
  image(z=out,col= grey (0:K/K),ylim=c(0,1),axes=F)
}

# Step 3. Visualize the results as an image
visu(H,3,0.1)
visu(H,4,0.1)
visu(H,5,0.1)
