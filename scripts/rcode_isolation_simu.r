# define the isolation function
isolation <- function(distance, area, alpha){
  n <- nrow(distance)
  isol <- rep(0,n)
  for(i in 1:n){
    isol[i] <- sum(area[-i] * exp(-alpha * distance[i,-i]))
  }
  return(isol)
}

# create some distance and area data
plots <- data.frame(x = runif(50),y = runif(50),frag_id = sample(1:50, 50, replace = FALSE))
area <- runif(50, 1, 100)
distt <- as.matrix(dist(plots[,c("x","y")]))

ii <- isolation(distt,area, 1) # compute the isolation value
isol <- ii[plots$frag_id]
isol <- scale(isol) # scale the isolation

X <- runif(50,-2,2) # a covariate
mu <- 1 + 2 * X - isol # the linear predictors
y <- rnorm(50,mu,0.1) # simulate some response
X <- matrix(c(rep(1,50),X),ncol=2,byrow=FALSE)
frag_id <- plots$frag_id

# fit the model
m <- stan(file = "~/Documents/PostDoc_Ghent/Spatialsynthesis_stuff/model/normal_isolation.stan",
          data = list(N=50,K=2,X=X,id_need=50,N_frag=50,y=y,frag_id=frag_id,distt=distt,area=area))
