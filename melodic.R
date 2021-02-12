###########################################################################
###########################################################################
# Model fitting: three different types
###########################################################################
###########################################################################


###########################################################################
# unconstrained: TYPE 1
###########################################################################

mldm1 = function(X, Y, M = 2, trace = FALSE, tol = 1e-8){
  # multivariate logistic distance model 
  # X: an I by P matrix with predictor variables 
  # Y: an I by R matrix with response variables
  # M: dimensionality of required solution
  # trace: if TRUE tracing information on the progress of the optimization is given
  # tol: convergence tolerance
  library(nnet)
  # preparation
  Yoriginal = Y
  Xoriginal = X
  I = nrow(Y)
  R = ncol(Y)
  P = ncol(X)
  # centering and scaling the predictor matrix. 
  X = scale(X)
  mx = attr(X, "scaled:center")
  sdx = attr(X, "scaled:scale")
  # create some structures
  # Y: super indicator matrix
  # Yhat: super matrix with estimated probabilities
  # J: block diagonal centering matrix
  # Ak/Al matrix for decomposition of V into K/Vk and L/Vl
  J.list <- YY <- vector("list", R)
  c.idx = c(0,rep(NA,R))
  Cr = rep(0,R) 
  for(r in 1:R){
    YY[[r]] = class.ind(Y[,r]) 
    c.idx[r+1] = c.idx[r] + ncol(YY[[r]])
    Cr[r] = ncol(YY[[r]])
    J.list[[r]] <- diag(Cr[r]) - 1/Cr[r] 
  }
  G = YY[[1]]; if(R>1){for(r in 2:R){G = cbind(G,YY[[r]])}}
  C = ncol(G)
  Ghat = matrix(NA, I, C)
  
  ones.I = matrix(1,I,1)
  
  J <- matrix(0, C, C)
  for (r in 1:R){
    j.ind <- (c.idx[r] + 1):c.idx[r + 1]
    J[j.ind, j.ind] <- J.list[[r]]
  }
  
  Ak = kronecker(diag(R), matrix(c(1,-1),2,1)) # for K - response variable discrimination
  Al = kronecker(diag(R), matrix(c(1,1),2,1)) # for L - response variable position
  
  # preparatory steps for rows
  eig.out = eigen(t(X) %*% X)
  Rx = eig.out$vectors %*% diag(sqrt(eig.out$values)) %*% t(eig.out$vectors)
  iXXX = solve(t(X) %*% X) %*% t(X)
  iRx = solve(Rx)
  iRxX = iRx %*% t(X)
    
  # initialization
  svd.out = mysvd2(Rx %*% iXXX %*% G, M, 0) 
  B = iRx %*% svd.out$u *sqrt(I)
  U = X %*% B
  V = svd.out$v %*% svd.out$D2 /sqrt(I)
  
  # compute initial estimates of K, L, Vk and Vl (in notitie patrick H, C, G1 en G2)
  K = (J %*% V)[seq(1,C,by = 2), ] # item discrimination
  L = ((diag(C) - J) %*% V)[seq(1,C,by = 2), ] # item locations
  Vk = Ak %*% K # deviation from midpoints/discrimination
  Vl = Al %*% L  # midpoints/location
  V = Vl + Vk
  
  # fitted values
  theta = - sqdist.r(U,V)
  for(r in 1:R){
    Ghat[,(c.idx[r]+1):c.idx[r+1]] = 
      exp(theta[,(c.idx[r]+1):c.idx[r+1]])/rowSums(exp(theta[,(c.idx[r]+1):c.idx[r+1]]))
  }
  
  # initial deviance
  QD = rep(NA,1e6)
  QD[1] = -2*sum(G * log(Ghat))
  
  #####################################################################
  #####################################################################
  # START ITERATING
  #####################################################################
  #####################################################################
  iter = 0; dif = 1
  while(dif > tol){
    iter = iter + 1
    # Update a
    # || .5 ZA + 1a' - XBK' ||^2
    Z = 1/2 * (theta + 2*(G - Ghat)) %*% Ak
    a = matrix(-colMeans(Z), R, 1)
    
    # if(M == 1){L = a/K}
    # if(M > 1){L = (a / rowSums(K)) %*% matrix(1,1,M)}
    # Vl = Al %*% L
    # V = Vk + Vl
    if(M == 1){L = a/K}
    if(M > 1){
      for(r in 1:R){
        L[r, ] = a[r] * K[r, ]/sum(K[r, ]^2)
      }
    }
    Vl = Al %*% L  # midpoints
    V = Vk + Vl
    
    
    # fitted values
    theta = - sqdist.r(U,V)
    for(r in 1:R){
      Ghat[,(c.idx[r]+1):c.idx[r+1]] = 
        exp(theta[,(c.idx[r]+1):c.idx[r+1]])/rowSums(exp(theta[,(c.idx[r]+1):c.idx[r+1]]))
    }
    
    # Update B and K
    Z = 1/2 * (theta + 2*(G - Ghat)) %*% Ak
    #Ztilde = iXXX %*% (Z + ones.I %*% t(a)) 
    #svd.out = mysvd2(Rx %*% Ztilde, M, 0)
    Ztilde = iRxX %*% (Z + ones.I %*% t(a)) 
    svd.out = mysvd2(Ztilde, M, 0)
    
    B = iRx %*% svd.out$u * sqrt(I)
    U = X %*% B
    K = svd.out$v %*% svd.out$D2 /sqrt(I)
    
    # Derive L, Vk and Vl (V = Vk + Vl)
    Vk = Ak %*% K # deviation from midpoints
    # if(M == 1){L = a/K}
    # if(M > 1){L = (a / rowSums(K)) %*% matrix(1,1,M)}
    # Vl = Al %*% L  # midpoints
    # V = Vk + Vl
    
    if(M == 1){L = a/K}
    if(M > 1){
      for(r in 1:R){
        L[r, ] = a[r] * K[r, ]/sum(K[r, ]^2)
      }
    }
    Vl = Al %*% L  # midpoints
    V = Vk + Vl
    
    # fitted values
    theta = - sqdist.r(U,V)
    for(r in 1:R){
      Ghat[,(c.idx[r]+1):c.idx[r+1]] = 
        exp(theta[,(c.idx[r]+1):c.idx[r+1]])/rowSums(exp(theta[,(c.idx[r]+1):c.idx[r+1]]))
    }
    
    # deviance
    QD[iter+1] = -2*sum(G * log(Ghat))
    dif = (QD[iter] - QD[iter+1])/QD[iter+1]
    if(trace){
      cat(iter, QD[iter], QD[iter+1], dif, "\n")
    }
  } # 
  #####################################################################
  #####################################################################
  # END ITERATING
  #####################################################################
  #####################################################################
  LRcoefs = matrix(NA, P, R)
  for(r in 1:R){
    dif = V[(r*2), , drop = FALSE] - V[((r-1)*2 +1), , drop = FALSE]
    LRcoefs[, r] = rowSums(B * matrix(1,P,1) %*% dif)
  }
  rownames(LRcoefs) = colnames(Xoriginal)
  colnames(LRcoefs) = colnames(Y)
  # finalisation
  results = list(
    Y = Y,
    Xoriginal = Xoriginal,
    G = G,
    X = X,
    mx = mx,
    sdx = sdx,
    B = B,
    U = U,
    V = V,
    a = a,
    K = K,
    L = L,
    Ghat = Ghat,
    QD = QD[(iter+1)],
    npar = (P + R)*M + R - M*(M+1)/2, 
    iter = iter, 
    J = J,
    LRcoef = LRcoefs
  )
  return(results)
}

###########################################################################
# constrained with unequal distances: TYPE 2
###########################################################################

mldm2 = function(X, Y, dim.indic, trace = FALSE, tol = 1e-8){
  # multivariate logistic distance model 
  # confirmatory weher we know which response variable belongs to which dimension
  # X: an I by P matrix with predictor variables 
  # Y: an I by R matrix with response variables
  # dim.indic: R x M matrix indicating which response goes to which dimension
  # trace: if TRUE tracing information on the progress of the optimization is given
  # tol: convergence tolerance
  library(nnet)
  # preparation
  Yoriginal = Y
  Xoriginal = X
  I = nrow(Y)
  R = ncol(Y)
  P = ncol(X)
  M = ncol(dim.indic)
  dim.indic = matrix(as.logical(dim.indic), R, M)
  # centering and scaling the predictor matrix. 
  X = scale(X)
  mx = attr(X, "scaled:center")
  sdx = attr(X, "scaled:scale")
  # create some structures
  # Y: super indicator matrix
  # Yhat: super matrix with estimated probabilities
  # J: block diagonal centering matrix
  # Ak/Al matrix for decomposition of V into K/Vk and L/Vl
  J.list <- YY <- vector("list", R)
  c.idx = c(0,rep(NA,R))
  Cr = rep(0,R) 
  for(r in 1:R){
    YY[[r]] = class.ind(Y[,r]) 
    c.idx[r+1] = c.idx[r] + ncol(YY[[r]])
    Cr[r] = ncol(YY[[r]])
    J.list[[r]] <- diag(Cr[r]) - 1/Cr[r] 
  }
  G = YY[[1]]; if(R>1){for(r in 2:R){G = cbind(G,YY[[r]])}}
  C = ncol(G)
  Ghat = matrix(NA, I, C)
  
  resp.indic = matrix(0, C, M)
  for(r in 1:R){
    resp.indic[(c.idx[r]+1): c.idx[r+1], ] = matrix(1, Cr[r], 1) %*% dim.indic[r, ]
  }
  resp.indic = matrix(as.logical(resp.indic), C, M)
  
  ones.I = matrix(1,I,1)
  
  J <- matrix(0, C, C)
  for (r in 1:R){
    j.ind <- (c.idx[r] + 1):c.idx[r + 1]
    J[j.ind, j.ind] <- J.list[[r]]
  }
  
  Ak = kronecker(diag(R), matrix(c(1,-1),2,1)) # for K - response variable discrimination
  Al = kronecker(diag(R), matrix(c(1,1),2,1)) # for L - response variable position
  
  # preparatory steps for rows
  eig.out = eigen(t(X) %*% X)
  Rx = eig.out$vectors %*% diag(sqrt(eig.out$values)) %*% t(eig.out$vectors)
  iXXX = solve(t(X) %*% X) %*% t(X)
  iRx = solve(Rx)
  iRxX = iRx %*% t(X)
  
  # initialization
  svd.out = mysvd2(iRxX %*% G, M, 0) 
  B = iRx %*% svd.out$u *sqrt(I)
  U = X %*% B
  V = resp.indic * svd.out$v %*% svd.out$D2 /sqrt(I)
  
  # compute initial estimates of K, L, Gk and Gl (in notitie patrick H, C, G1 en G2)
  K = (J %*% V)[seq(1,C,by = 2), ] # item discrimination
  L = ((diag(C) - J) %*% V)[seq(1,C,by = 2), ] # item locations
  K = dim.indic * ifelse(K > 0.05, K, 0.05)
  
  Vk = Ak %*% K # deviation from midpoints/discrimination
  Vl = Al %*% L  # midpoints/location
  V = Vl + Vk
  
  # fitted values
  theta = - sqdist.r(U,V)
  for(r in 1:R){
    Ghat[,(c.idx[r]+1):c.idx[r+1]] = 
      exp(theta[,(c.idx[r]+1):c.idx[r+1]])/rowSums(exp(theta[,(c.idx[r]+1):c.idx[r+1]]))
  }
  
  # initial deviance
  QD = rep(NA,1e6)
  QD[1] = -2*sum(G * log(Ghat))
  
  #####################################################################
  #####################################################################
  # START ITERATING
  #####################################################################
  #####################################################################
  iter = 0; dif = 1
  while(dif > tol){
    iter = iter + 1
    # Update a
    Z = 1/2 * (theta + 2*(G - Ghat)) %*% Ak
    a = matrix(-colMeans(Z), R, 1)
    L = dim.indic * ((a / rowSums(dim.indic * K)) %*% matrix(1,1,M))
    
    Vl = Al %*% L
    V = Vk + Vl
    
    # fitted values
    theta = - sqdist.r(U,V)
    for(r in 1:R){
      Ghat[,(c.idx[r]+1):c.idx[r+1]] = 
        exp(theta[,(c.idx[r]+1):c.idx[r+1]])/rowSums(exp(theta[,(c.idx[r]+1):c.idx[r+1]]))
    }
    
    # Update B and K - dimensionwise
    for(m in 1:M){
      Z = 1/2 * (theta + 2*(G - Ghat)) %*% Ak
      Ztilde = iRxX %*% (Z + ones.I %*% t(a) - U[, -m] %*% t(K[, -m])) 
      svd.out = mysvd2(Ztilde[, dim.indic[ ,m]], 1, 0)
      B[ ,m] = iRx %*% svd.out$u * sqrt(I)
      U = X %*% B
      K[dim.indic[, m], m] = svd.out$v %*% svd.out$D2 /sqrt(I)
      
      # Derive L, Vk and Vl (V = Vk + Vl)
      Vk = Ak %*% K # deviation from midpoints
      L = dim.indic * ((a / rowSums(dim.indic * K)) %*% matrix(1,1,M))
      Vl = Al %*% L  # midpoints
      V = Vk + Vl
      
      # fitted values
      theta = - sqdist.r(U,V)
      for(r in 1:R){
        Ghat[,(c.idx[r]+1):c.idx[r+1]] = 
          exp(theta[,(c.idx[r]+1):c.idx[r+1]])/rowSums(exp(theta[,(c.idx[r]+1):c.idx[r+1]]))
      }
    }
    # deviance
    QD[iter+1] = -2*sum(G * log(Ghat))
    dif = (QD[iter] - QD[iter+1])/QD[iter+1]
    if(trace){
      cat(iter, QD[iter], QD[iter+1], dif, "\n")
    }
  } # 
  #####################################################################
  #####################################################################
  # END ITERATING
  #####################################################################
  #####################################################################
  LRcoefs = matrix(NA, P, R)
  for(r in 1:R){
    dif = V[(r*2), , drop = FALSE] - V[((r-1)*2 +1), , drop = FALSE]
    LRcoefs[, r] = rowSums(B * matrix(1,P,1) %*% dif)
  }
  rownames(LRcoefs) = colnames(Xoriginal)
  colnames(LRcoefs) = colnames(Y)
  # finalisation
  results = list(
    Y = Y,
    Xoriginal = Xoriginal,
    G = G,
    X = X,
    mx = mx,
    sdx = sdx,
    B = B,
    U = U,
    V = V,
    a = a,
    K = K,
    L = L,
    Ghat = Ghat,
    QD = QD[(iter+1)],
    npar = (P-1)*M + R + sum(dim.indic),
    iter = iter, 
    J = J,
    LRcoef = LRcoefs
  )
  return(results)
}

###########################################################################
# constrained with equal distances: TYPE 3
###########################################################################

mldm3 = function(X,Y,dim.indic){
  # 
  I = nrow(X)
  P = ncol(X)
  R = ncol(Y)
  M = ncol(dim.indic)
  
  Xoriginal = X
  Yoriginal = Y
  
  # standardizing X
  XX = scale(X)
  mx = attr(XX, "scaled:center")
  sx = attr(XX, "scaled:scale")
  X = matrix(XX, I, P)
  
  # Y indicator
  YY <- vector("list", R)
  c.idx = c(0,rep(NA,R))
  Cr = rep(0,R) 
  for(r in 1:R){
    YY[[r]] = class.ind(Y[,r]) 
    c.idx[r+1] = c.idx[r] + ncol(YY[[r]])
    Cr[r] = ncol(YY[[r]])
  }
  G = YY[[1]]; if(R>1){for(r in 2:R){G = cbind(Y,YY[[r]])}}
  Ghat = matrix(NA, nrow(Y), ncol(G))
  C = ncol(G)
  
  # transforming to glm-structure
  yvec = matrix(t(Yoriginal), I*R, 1)
  XX = matrix(NA, I*R, (R + P*M))
  for(i in 1:I){
    XX[((i-1)*R + 1):(i*R),] = cbind(diag(R), dim.indic %x% Xoriginal[i,,drop = FALSE])
  }
  #
  dim.sum = rowSums((dim.indic))
  glm.output = glm.fit(XX, yvec, family = binomial(), intercept = FALSE)
  V = matrix(0,(2*R),M)
  for (r in 1:R){
    for (m in 1:M){
      if(dim.indic[r,m] == 1){
        # check if loads on more dims - dim.sum
        V[((r-1)*2 + 1),m] = - (glm.output$coef[r] + 0.5)
        V[(r*2),m] = V[((r-1)*2 + 1),m] + 1
      }
    }
  }
  B = matrix(glm.output$coef[(R+1):ncol(XX)],P,M, byrow = FALSE)
  U = X %*% B
  theta = -sqdist.r(U,V)
  Ghat = matrix(NA, nrow(G), ncol(G))
  for(r in 1:R){
    Ghat[,(c.idx[r]+1):c.idx[r+1]] =   exp(theta[,(c.idx[r]+1):c.idx[r+1]])/rowSums(exp(theta[,(c.idx[r]+1):c.idx[r+1]]))
  }
  LRcoefs = matrix(NA, P, R)
  for(r in 1:R){
    dif = V[(r*2), , drop = FALSE] - V[((r-1)*2 +1), , drop = FALSE]
    LRcoefs[, r] = rowSums(B * matrix(1,P,1) %*% dif)
  }
  rownames(LRcoefs) = colnames(Xoriginal)
  colnames(LRcoefs) = colnames(Y)
  
  # finalisation
  results = list(
    Y = Y,
    Xorginal = X,
    G = G,
    X = X,
    mx = mx,
    sdx = sdx,
    B = B,
    U = U,
    V = V,
    a = NULL,
    K = NULL,
    L = NULL,
    Ghat = Ghat,
    QD = glm.output$deviance,
    npar = length(glm.output$coef),
    iter = glm.output$iter,
    J = NULL,
    LRcoef = LRcoefs
  )
  return(results)  
}



###########################################################################
###########################################################################
# visualisation
###########################################################################
###########################################################################

plot.mldm = function(output, var.axis = TRUE, dec.lines = TRUE){
  # plotting function for 2D outcome object
  library(ggplot2)
  library(ggforce)
  library(ggrepel)
  library(directlabels)
  
  xnames = colnames(output$X)
  ynames = colnames(output$Y)
  ynames = paste0(rep(ynames, each = 2), c(0,1))
  
  U = as.data.frame(output$U)
  colnames(U) = c("Dim1", "Dim2")
  
  V = as.data.frame(output$V)
  colnames(V) = c("Dim1", "Dim2")
  rownames(V) = ynames
  
  X = output$X
  P = ncol(X)
  B = output$B
  R = ncol(output$Y)
  
  maxX = apply(output$X,2,max)
  BB = maxX*B
  BB = as.data.frame(BB)
  colnames(BB) = c("Dim1", "Dim2")
  rownames(BB) = xnames
  
  # getting labels outside borders
  idx1 = apply(abs(B), 1, which.max)
  t = s = rep(NA,P)
  for(pp in 1:P){
    t[pp] = 3.4/(abs(B[pp,idx1[pp]])) * B[pp,-idx1[pp]]
    s[pp] = sign(B[pp,idx1[pp]])
  }
  CC = cbind(idx1, t, s)
  
  # plot category points and person points
  p = ggplot() +
    geom_point(data = V, aes(x = Dim1, y = Dim2), colour = "darkgreen", size = 3) +
    geom_point(data = U, aes(x = Dim1, y = Dim2), shape = ".") + 
    xlab("Dimension 1") + 
    ylab("Dimension 2") + 
    coord_fixed()   
    #theme(axis.text.x = element_blank(),axis.text.y = element_blank())
  # add variable axes
  if(var.axis){
    p = p + geom_abline(intercept = 0, slope = B[,2]/B[,1], colour = "lightskyblue")
    # add markers to variable axis - gebruik CC
    # if (idx1 = 1) & s = 1 => label secondary y axis
    # if (idx1 = 1) & s = -1 => label primary y axis
    # if (idx1 = 2) & s= 1 => label secondary x axis
    # if (idx1 = 2) & s= -1 => label primary x axis
    bottom = which(CC[, "idx1"] == 2 & CC[, "s"] == -1)
    top =  which(CC[, "idx1"] == 2 & CC[, "s"] == 1)
    right = which(CC[, "idx1"] == 1 & CC[, "s"] == 1)
    left = which(CC[, "idx1"] == 1 & CC[, "s"] == -1)
    
    p = p + scale_x_continuous(limits = c(-3,3), breaks = CC[bottom, "t"], labels = xnames[bottom], 
                               sec.axis = sec_axis(trans ~ ., breaks = CC[top, "t"], labels = xnames[top]))
    p = p + scale_y_continuous(limits = c(-3,3), breaks = CC[left, "t"], labels = xnames[left], 
                               sec.axis = sec_axis(trans ~ ., breaks = CC[right, "t"], labels = xnames[right]))
    # oude variable markers
    # p = p + annotate("text", x= BB[,1], y = BB[,2], label = xnames, angle = 180/pi*atan(B[,2]/B[,1]), hjust = sign(-BB[,1]))
    # add variable markers
    for (r in 1:P){
      B1 = data.frame(outer(seq(ceiling(min(X[,r])),floor(max(X[,r]))),B[r,])) 
      colnames(B1) = c("Dim1", "Dim2")
      rownames(B1) = as.character(seq(ceiling(min(X[,r])),floor(max(X[,r]))))
      p = p + geom_point(data = B1, aes(x= Dim1, y= Dim2), colour = "lightskyblue", size = 2)
      #p = p + geom_text(data = B1, aes(x= Dim1, y= Dim2, label = rownames(B1), family = "mono"))
    }
  }
  
  
  # decision lines
  if(dec.lines){
    int.slp = matrix(NA,R,2)
    mid = matrix(NA,R,2)
    my.slp = rep(NA, R)
    VV = as.matrix(V)
    for (r in 1:R){
      mid[r, ] = 0.5 * (VV[(r*2), ] + VV[(r*2-1), ])
      dif = V[(r*2), , drop = FALSE]-V[(r*2-1), , drop = FALSE]
      my.slp[r] = -(dif[1, 1]/dif[1, 2])   
      my.int = mid[r, 2] - mid[r, 1] * my.slp[r]
      int.slp[r,1] = my.int
      int.slp[r,2] = my.slp[r]
    }
    
    int.slp = as.data.frame(int.slp)
    idx = is.finite(my.slp)
    p = p + geom_abline(intercept = int.slp[idx,1], slope = int.slp[idx,2], colour = "darkgreen") 
    p = p + geom_vline(xintercept = mid[!idx,1], colour = "darkgreen")
  }
  # 
  p = p + geom_text_repel(data = V, aes(x= Dim1, y = Dim2, label = ynames, family = "mono")) 
  p = p + theme(text=element_text(size=12,  family="mono"))
  p = p + theme_classic()
  return(p)
}

mldm.diag = function(output){
  # diagnosis of mldm 
  # input: mldm output object
  # return: Quality of Representation
  Y = output$Y
  X = output$Xoriginal
  R = ncol(Y)
  devs.lr = devs.null = rep(NA,R)
  for(r in 1:R){
    v = glm(Y[,r] ~ X, family = "binomial")
    devs.null[r] = v$null.deviance
    devs.lr[r] = v$deviance
  }
  devs.mldm = rowSums(matrix(-2*colSums(output$G * log(output$Ghat)), R, 2, byrow = TRUE)) 
  # quality of representation
  QOR = (devs.null - devs.mldm) /(devs.null - devs.lr)
  QoR = matrix(QOR, R, 1)
  rownames(QoR) = colnames(Y)
  return(QoR)
}

###########################################################################
###########################################################################
# cross validation
###########################################################################
###########################################################################

predict.mldm = function(output, newX, newY = NULL){
  # predict function for mldm models
  # outout is an output object from mldm1, mldm2, or mldm3
  # newX is a matrix with new predictor variable values
  # newY is optional matrix with responses for cases in newX
  
  X = scale(newX, center = output$mx, scale = output$sdx)
  B = output$B
  V = output$V
  U = X %*% B
  
  R = ncol(output$Y)
  n = nrow(newX)
  
  Ghat = matrix(NA, n, 2*R)
  if(!is.null(newY)){
    YY <- vector("list", R)
    for(r in 1:R){
      YY[[r]] = cbind(1 - newY[ , r], newY[ , r]) 
    }
    G = YY[[1]]; if(R>1){for(r in 2:R){G = cbind(G,YY[[r]])}} 
  }

  for(r in 1:R){
    Vr = V[(2*r-1):(2*r),]
    D2 = (outer(diag(U %*% t(U)), rep(1,2)) + outer(rep(1,n), diag(Vr %*%t(Vr))) - 2 * U %*% t(Vr))/2
    Ghat[,(2*r-1):(2*r)] = cbind(exp(-D2[,1])/(exp(-D2[,1])+exp(-D2[,2])), exp(-D2[,2])/(exp(-D2[,1])+exp(-D2[,2])))
  }
  
  # performance measures
  # xval deviance
  # xval brier score
  if(!is.null(newY)){
    
    # make two list for response specific performances
    dev.xval.r <- vector(mode = "list", length = length(R)) 
    brier.xval.r <- vector(mode = "list", length = length(R)) 
    # fill lists
    for(r in 1:R){
      dev.xval.r[[r]] = -2*sum(G[ , (2*r-1):(2*r)] * log(Ghat[ , (2*r-1):(2*r)]))
      brier.xval.r[[r]] = sqrt(sum(((G - Ghat)[ , 2*r])^2)/n)  
    }

    # overall performance
    dev.xval = Reduce( "+", dev.xval.r) # total deviance 
    brier.xval = Reduce( "+" , brier.xval.r)/R # average brier score
  }
  
  # output object
  if(is.null(newY)){
    result = list(
    output = output,
    newX = newX,
    X = X,
    newY = NULL,
    Ghat = Ghat,
    G = NULL,
    dev.xval = NULL,
    brier.xval = NULL,
    dev.xval.r = NULL,
    brier.xval.r = NULL   
    )
  }
  else {
    result = list(
      output = output,
      newX = newX,
      X = X,
      newY = newY,
      Ghat = Ghat,
      G = G,
      dev.xval = dev.xval,
      brier.xval = brier.xval,
      dev.xval.r = dev.xval.r,
      brier.xval.r = brier.xval.r   
    )
  }
  return(result)
}



###########################################################################
###########################################################################
# help functions
###########################################################################
###########################################################################

sqdist.r = function(U,V){
  # computes the matrix with ``reduced'' squared Euclidean distances
  ones.I = matrix(1,nrow(U),1)
  D =  ones.I %*% t(diag(V %*% t(V))) - 2* U %*% t(V)
  return(D/2)
}

mysvd = function(A, M){
  # an SVD procedure that keeps matrices
  I = nrow(A)
  J = ncol(A)
  svd.out = svd(A, nu = M, nv = M)
  U = matrix(svd.out$u, I, M)
  V = matrix(svd.out$v, J, M)
  D = diag(svd.out$d, M, M)
  result = list(
    u = U,
    v = V,
    d = D
  )
  return(result)
}

mysvd2 = function(A, M, tau){
  # an SVD procedure that keeps matrices
  I = nrow(A)
  J = ncol(A)
  svd.out = svd(A, nu = M, nv = M)
  U = matrix(svd.out$u, I, M)
  V = matrix(svd.out$v, J, M)
  D = diag(svd.out$d, M, M)
  D1 = diag(svd.out$d^tau, M, M)
  D2 = diag(svd.out$d^(1-tau), M, M)
  result = list(
    u = U,
    v = V,
    d = D,
    D1 = D1,
    D2 = D2
  )
  return(result)
}

