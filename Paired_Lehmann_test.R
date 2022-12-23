# Test the hypothesis that two dependent variables (x,y) have equal variances.
# Paired Lehmann-type statistic estimate the probability P(|x1 - x2| < |y1 - y2|).

paired.Lehmann.Var.test <- function(x, y = NULL, method = c("NAP", "PB"), alternative = c("two.sided","less","greater"), nboot = 1000){
  ### Paired Lehmann-type test to assess the difference between two scale parameters in paired data
  # x is first sample.
  # y is second sample.
  # "NAP" denotes the p-value based on the normal approximation (NAP)
  # "PB" denotes the p-value based on the percentile bootstrap (PB)
  # nboot is the number of bootstrap samples. 
  
  alternative <- match.arg(alternative)  ## alternative hypothesis
  method <- match.arg(method)            ## approximation method of p-value
  if(!is.null(y)){
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  }else{
    stop("'y' is missing for paired test")
  }

  n <- length(x)            ## Sample size
  X <- abs(outer(x,x,"-"))  ## First sample 
  Y <- abs(outer(y,y,"-"))  ## Second sample
  sign.diff <- (sign(Y - X) + 1)/2              
  L <- mean(sign.diff[upper.tri(sign.diff)])
  
  if(method == "NA"){
    method.name <- c("Paired Lehmann-type test based on the normal approximation")
    SUM <- function(i) (sum(sign.diff[i,]) - 0.5)^2 
    L.Var.est <- 4*(n-2)/(n*(n-1))*(sum(sapply(1:n,SUM))/(n*(n-1)^2) - L^2) + 1/(2*n*(n-1))
    if(alternative == "two.sided"){
      pval <- 2 * pnorm(abs((L - 0.5)/sqrt(L.Var.est)), lower.tail = FALSE)
    }
    if(alternative == "less"){
      pval <- pnorm((L - 0.5)/sqrt(L.Var.est), lower.tail = FALSE)
    }
    if(alternative == "greater"){
      pval <- pnorm((L - 0.5)/sqrt(L.Var.est), lower.tail = TRUE)
    }
  }
  if(method == "PB"){
    method.name <- c("Paired Lehmann-type test based on the percentile bootstrap")
    stat.boot <- numeric(nboot)
    z <- cbind(x,y)
    for(b in 1:nboot){
      zb <- z[sample(1:n, replace = TRUE),]
      Xb <- abs(outer(zb[,1],zb[,1],"-")) 
      Yb <- abs(outer(zb[,2],zb[,2],"-")) 
      sign.diff.b <- (sign(Yb - Xb) + 1)/2              
      stat.boot[b] <- mean(sign.diff.b[upper.tri(sign.diff.b)])
    }
    if(alternative == "two.sided"){
      pval <- 2 * min(sum(stat.boot - 0.5 < 0)/nboot, sum(stat.boot - 0.5 > 0)/nboot)
    }
    if(alternative == "less"){
      pval <- sum(stat.boot - 0.5 < 0)/nboot
    }
    if(alternative == "greater"){
      pval <- sum(stat.boot - 0.5 > 0)/nboot
    }
  }
  names(L) <- "Lehmann.stat"
  rval <- list(statistic = L, p.value = pval, alternative = alternative, method = method.name, data.name = dname)
  class(rval) <- "htest"
  return(rval)
}

## Example 1  (Bick R, Adams T, Schmalhorst W. (1976). Bleeding times, platelet adhesion, and aspirin. American Journal of Clinical Pathology. 65: 69–72)
X <- matrix(c(270, 150, 270, 420, 202, 255, 165, 220, 305, 210, 240, 300, 300,  70, 265, 215,  95,  85, 200,
              525, 570, 190, 395, 370, 210, 490, 250, 360, 285, 630, 385, 195, 295, 275, 270, 295, 535, 345), nrow = 19, ncol = 2)
paired.Lehmann.Var.test(X[,1],X[,2], method = "NAP", alternative = "less")                  # p-value = 0.05858
paired.Lehmann.Var.test(X[,1],X[,2], method = "PB",  alternative = "less",nboot = 1000000)  # p-value = 0.06187 if set.seed(1234)


## Example 2 (Bick R, Adams T, Schmalhorst W. (1976). Bleeding times, platelet adhesion, and aspirin. American Journal of Clinical Pathology. 65: 69–72)
X <- matrix(c(67, 74, 90, 82, 90, 82, 62, 76, 89, 84, 83, 94, 84, 69, 76, 77, 72, 93, 89,
              67, 82, 73, 80, 78, 84, 57, 63, 93, 89, 83, 97, 90, 40, 89, 61, 90, 92, 84), nrow = 19, ncol = 2)
paired.Lehmann.Var.test(X[,1],X[,2], method = "NAP", alternative = "less")                  # p-value = 0.04886
paired.Lehmann.Var.test(X[,1],X[,2], method = "PB",  alternative = "less",nboot = 1000000)  # p-value = 0.05183 if set.seed(1234)
