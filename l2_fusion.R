# L2 Fusion approach optimization


#' Generate block diagonal matrices to allow for fused L2 optimization with glmnet.
#'
#' @param X covariates matrix (n by p).
#' @param Y response vector (length n).
#' @param groups vector of group indicators (ideally factors, length n)
#' @param G matrix representing the fusion strengths between pairs of
#' groups (K by K). Zero entries are assumed to be independent pairs.
#' @param intercept whether to include an (per-group) intercept in the model
#' @param penalty.factors vector of weights for the penalization of
#' each covariate (length p)
#' @param scaling Whether to scale each subgroup by its size. See Details for an explanation.
#'
#'
#' @return A list with components X, Y, X.fused and penalty, where
#' X is a n by pK block-diagonal bigmatrix, Y is a
#' re-arranged bigvector of length n, and X.fused is a
#' choose(K,2)*p by pK bigmatrix encoding the fusion penalties.
#'
#' @details We use the \code{glmnet} package to perform fused subgroup regression.
#' In order to achieve this, we need to reformulate the problem as Y' = X'beta',
#' where Y' is a concatenation of the responses Y and a vector of zeros, X' is a
#' a matrix consisting of the block-diagonal matrix n by pK matrix X, where each
#' block contains the covariates for one subgroups, and the choose(K,2)*p by pK
#' matrix encoding the fusion penalties between pairs of groups. The vector of
#' parameters beta' of length pK can be rearranged as a p by K matrix giving the
#' parameters for each subgroup. The lasso penalty on the parameters is handled
#' by glmnet.
#'
#' One weakness of the approach described above is that larger subgroups will
#' have a larger influence on the global parameters lambda and gamma.
#' In order to mitigate this, we introduce the \code{scaling} parameter. If
#' \code{scaling=TRUE}, then we scale the responses and covariates for each
#' subgroup by the number of samples in that group.
#'
#'
#' @import Matrix
#'
#' @export
#'
#' @examples
#' set.seed(123)
#'
#' # Generate simple heterogeneous dataset
#' k = 4 # number of groups
#' p = 100 # number of covariates
#' n.group = 15 # number of samples per group
#' sigma = 0.05 # observation noise sd
#' groups = rep(1:k, each=n.group) # group indicators
#' # sparse linear coefficients
#' beta = matrix(0, p, k)
#' nonzero.ind = rbinom(p*k, 1, 0.025/k) # Independent coefficients
#' nonzero.shared = rbinom(p, 1, 0.025) # shared coefficients
#' beta[which(nonzero.ind==1)] = rnorm(sum(nonzero.ind), 1, 0.25)
#' beta[which(nonzero.shared==1),] = rnorm(sum(nonzero.shared), -1, 0.25)
#'
#' X = lapply(1:k,
#'            function(k.i) matrix(rnorm(n.group*p),
#'                                 n.group, p)) # covariates
#' y = sapply(1:k,
#'            function(k.i) X[[k.i]] %*% beta[,k.i] +
#'                            rnorm(n.group, 0, sigma)) # response
#' X = do.call('rbind', X)
#'
#' # Pairwise Fusion strength hyperparameters (tau(k,k'))
#' # Same for all pairs in this example
#' G = matrix(1, k, k)
#'
#' # Generate block diagonal matrices
#' transformed.data = generateBlockDiagonalMatrices(X, y, groups, G)
generateBlockDiagonalMatrices <- function(X, Y, groups, G, intercept=FALSE,
                                          penalty.factors=rep(1, dim(X)[2]),
                                          scaling=FALSE) {
  group.names = sort(unique(groups))
  num.groups = length(group.names)
  
  num.pairs = num.groups*(num.groups-1)/2
  
  if(intercept) X = cbind(X, matrix(1, dim(X)[1], 1)) # Include intercept
  
  new.y = Matrix(0, length(Y)+num.pairs*dim(X)[2],1, sparse=TRUE)
  new.x = Matrix(0, dim(X)[1], dim(X)[2]*num.groups,
                 sparse=TRUE)
  new.x.f = Matrix(0, num.pairs*dim(X)[2], dim(X)[2]*num.groups,
                   sparse=TRUE)
  
  row.start = 1
  col.start = 1
  
  # Add data into block matrix
  for(group.i in 1:num.groups) {
    group.inds = groups==group.names[group.i]
    
    row.range = row.start:(row.start+sum(group.inds)-1)
    col.range = col.start:(col.start+dim(X)[2]-1)
    
    new.y[row.range] = Y[group.inds]
    new.x[row.range, col.range] = X[group.inds,]
    
    if(scaling) {
      scale_factor = sqrt(sum(group.inds))
      new.y[row.range] = new.y[row.range]/scale_factor
      new.x[row.range, col.range] = new.x[row.range, col.range]/scale_factor
    }
    
    row.start = row.start + sum(group.inds)
    col.start = col.start + dim(X)[2]
  }
  
  col.start.i = 1
  row.start = 1
  
  # Add gamma contraints into block matrix
  for(group.i in 1:(num.groups-1)) {
    
    col.start.j = col.start.i + dim(X)[2]
    
    col.range.i = col.start.i:(col.start.i+dim(X)[2]-1)
    
    for(group.j in (group.i+1):num.groups) {
      
      tau = G[group.i, group.j]
      
      row.range = row.start:(row.start+dim(X)[2]-1)
      col.range.j = col.start.j:(col.start.j+dim(X)[2]-1)
      
      new.x.f[cbind(row.range, col.range.i)] = sqrt(tau)
      new.x.f[cbind(row.range, col.range.j)] = -sqrt(tau)
      
      if(intercept) { # Don't fuse intercept
        new.x.f[row.range[length(row.range)],
                col.range.i[length(col.range.i)]] = 0
        new.x.f[row.range[length(row.range)],
                col.range.j[length(col.range.j)]] = 0
      }
      
      row.start = row.start + dim(X)[2]
      col.start.j = col.start.j + dim(X)[2]
    }
    
    col.start.i = col.start.i + dim(X)[2]
  }
  
  if(intercept) {
    penalty = rep(c(penalty.factors, 0), num.groups)
  } else {
    penalty = rep(penalty.factors, num.groups)
  }
  
  return(list(X=new.x, Y=new.y, X.fused=new.x.f, penalty=penalty))
}

#' Optimise the fused L2 model with glmnet (using transformed input data)
#'
#' @param transformed.x Transformed covariates (output of generateBlockDiagonalMatrices)
#' @param transformed.x.f Transformed fusion constraints (output of generateBlockDiagonalMatrices)
#' @param transformed.y Transformed response (output of generateBlockDiagonalMatrices)
#' @param groups Grouping factors for samples (a vector of size n, with K factor levels)
#' @param lambda Sparsity penalty hyperparameter
#' @param gamma Fusion penalty hyperparameter
#' @param ... Further options passed to glmnet.
#'
#' @return Matrix of fitted beta values.
#' @import glmnet
#' @export
#' @return
#' A matrix with the linear coefficients for each group (p by k).
#'
#' @examples
#'
#' #' set.seed(123)
#'
#' # Generate simple heterogeneous dataset
#' k = 4 # number of groups
#' p = 100 # number of covariates
#' n.group = 15 # number of samples per group
#' sigma = 0.05 # observation noise sd
#' groups = rep(1:k, each=n.group) # group indicators
#' # sparse linear coefficients
#' beta = matrix(0, p, k)
#' nonzero.ind = rbinom(p*k, 1, 0.025/k) # Independent coefficients
#' nonzero.shared = rbinom(p, 1, 0.025) # shared coefficients
#' beta[which(nonzero.ind==1)] = rnorm(sum(nonzero.ind), 1, 0.25)
#' beta[which(nonzero.shared==1),] = rnorm(sum(nonzero.shared), -1, 0.25)
#'
#' X = lapply(1:k,
#'            function(k.i) matrix(rnorm(n.group*p),
#'                                 n.group, p)) # covariates
#' y = sapply(1:k,
#'            function(k.i) X[[k.i]] %*% beta[,k.i] +
#'                            rnorm(n.group, 0, sigma)) # response
#' X = do.call('rbind', X)
#'
#' # Pairwise Fusion strength hyperparameters (tau(k,k'))
#' # Same for all pairs in this example
#' G = matrix(1, k, k)
#'
#' # Generate block diagonal matrices
#' transformed.data = generateBlockDiagonalMatrices(X, y, groups, G)
#'
#' # Use L2 fusion to estimate betas (with near-optimal information
#' # sharing among groups)
#' beta.estimate = fusedL2DescentGLMNet(transformed.data$X,
#'                                      transformed.data$X.fused,
#'                                      transformed.data$Y, groups,
#'                                      lambda=c(0,0.001,0.1,1),
#'                                      gamma=0.001)
fusedL2DescentGLMNet <- function(X, y, groups, lambda = NULL, G = NULL, gamma=1, scaling = FALSE,...) {
  #use_scaling = (length(unique(table(groups))) > 1)
  
  # Smart default for G
  if(is.null(G)) {
    k = length(sort(unique(groups)))
    G = matrix(1, k, k)
  }
  
  transformed.data = generateBlockDiagonalMatrices(X, y, groups, G, scaling = scaling)
  transformed.x = transformed.data$X
  transformed.x.f = transformed.data$X.fused
  transformed.y = transformed.data$Y
  
  # Incorporate fusion penalty global hyperparameter
  transformed.x.f = transformed.x.f * sqrt(gamma * (dim(transformed.x)[1] + dim(transformed.x.f)[1]))
  transformed.x = rbind(transformed.x, transformed.x.f)
  
  group.names = sort(unique(groups))
  num.groups = length(group.names)
  
  # Change: No loop needed since lambda is a single value
  # Adjust correction factor based on scaling
  if(scaling) {
    correction_factor = length(unique(groups)) / dim(transformed.x)[1]
  } else {
    correction_factor = length(groups) / dim(transformed.x)[1]
  }
  
  # Replace glmnet with cv.glmnet using minimal folds
  #cv.result = cv.glmnet(transformed.x, transformed.y, 
      #                  standardize=FALSE,
     #                   nlambda = 10,
    #                    lambda_min_ratio = 0.001,
   #                     nfolds=3, ...)
  #corrected_lambda = cv.result$lambda.1se * correction_factor
  #coef.temp = coef(cv.result, s=corrected_lambda)
  # Keep original matrix structure
  #beta.mat = matrix(coef.temp[2:length(coef.temp)], 
    #                dim(transformed.x)[2]/num.groups, 
   #                 num.groups)
  #return(beta.mat)
  
  glmnet.result = glmnet(transformed.x, 
                         transformed.y, 
                         standardize=FALSE,
                         intercept = FALSE,
                         lambda=lambda*correction_factor, ...)
  beta.mat = matrix(NA, dim(transformed.x)[2]/num.groups, num.groups)
  
  coef.temp = coef(glmnet.result,
                   s=lambda*correction_factor ) # Correction for extra dimensions
  beta.mat = matrix(coef.temp[2:length(coef.temp)], 
                   dim(transformed.x)[2]/num.groups, 
                    num.groups)
  
  return(beta.mat)
}



#' Predict using a fitted fused L2 model
#'
#' @param beta.mat Matrix of fitted coefficients (p by k) from fusedL2DescentGLMNet
#' @param newX New covariates matrix (n_new by p)
#' @param newGroups Group indicators for new data (length n_new)
#' @param groups Original group indicators used in training
#' @param intercept Whether the model includes an intercept
#'
#' @return Predicted values (length n_new)
#' @export
#'
#' @examples
#' # After fitting the model
#' beta.estimate = fusedL2DescentGLMNet(X, y, groups, G, lambda=0.1, gamma=0.01)
#' 
#' # Make predictions on new data
#' predictions = predictFusedL2(beta.estimate, X_new, groups_new, groups)
predictFusedL2 <- function(beta.mat, newX, newGroups, groups, intercept=FALSE) {
  
  # Get unique group names from training
  train.group.names = sort(unique(groups))
  new.group.names = sort(unique(newGroups))
  
  # Check if all new groups were seen in training
  unseen.groups = setdiff(new.group.names, train.group.names)
  if(length(unseen.groups) > 0) {
    warning(paste("Unseen groups in new data:", paste(unseen.groups, collapse=", "),
                  "\nUsing average of all groups for these predictions"))
  }
  
  # Initialize predictions
  n.new = nrow(newX)
  predictions = numeric(n.new)
  
  # If intercept was used, need to handle it
  if(intercept) {
    # The last row of beta.mat contains intercepts
    p = nrow(beta.mat) - 1
    intercepts = beta.mat[p+1, ]
    beta.mat = beta.mat[1:p, , drop=FALSE]
    newX = cbind(newX, 1)  # Add intercept column
  }
  
  # Make predictions for each group
  for(group.i in new.group.names) {
    group.idx = which(newGroups == group.i)
    
    if(group.i %in% train.group.names) {
      # Use the specific group's coefficients
      group.col = which(train.group.names == group.i)
      
      if(intercept) {
        predictions[group.idx] = newX[group.idx, 1:p] %*% beta.mat[, group.col] + 
          intercepts[group.col]
      } else {
        predictions[group.idx] = newX[group.idx, ] %*% beta.mat[, group.col]
      }
    } else {
      # For unseen groups, use average of all groups
      avg.beta = rowMeans(beta.mat)
      
      if(intercept) {
        avg.intercept = mean(intercepts)
        predictions[group.idx] = newX[group.idx, 1:p] %*% avg.beta + avg.intercept
      } else {
        predictions[group.idx] = newX[group.idx, ] %*% avg.beta
      }
    }
  }
  
  return(predictions)
}
