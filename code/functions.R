library(Rcpp)
sourceCpp("/blue/amolstad/jinwen.fu/effect_aggregation/reproduce/code/algos.cpp")
library(dplyr)
library(igraph)
library(MASS)
library(rare)
library(fossil)
library(mclust)
library(profvis)
library(glmnet)
library(pROC)


############################## This section handles the trees #################################

# Calculates the depth (distance from root) of a given node in a parent-referenced tree.
# Input: node (scalar node id), df (data.frame with at least columns: node, parent).
# Output: integer depth (0 for root, 1 for children of root, etc.).
calculate_depth <- function(node, df) {
  depth <- 0
  while (!is.na(node)) {
    node <- df$parent[df$node == node]
    depth <- depth + 1
  }
  return(depth - 1)
}

# Identifies all non-leaf nodes and computes their depths in the tree.
# Input: df (data.frame with at least columns: node, parent).
# Output: data.frame with columns: node (non-leaf ids), depth (integer depths).
detect_non_leaf_nodes <- function(df) {
  non_leaf_nodes <- unique(df$parent[!is.na(df$parent)])
  non_leaf_depths <- sapply(non_leaf_nodes, calculate_depth, df = df)
  non_leaf_df <- data.frame(node = non_leaf_nodes, depth = non_leaf_depths)
  
  return(non_leaf_df)
}

# Collects all descendant leaves under a given node using a BFS over parent links.
# Input: node (scalar node id), df (data.frame with at least columns: node, parent).
# Output: sorted integer vector of leaf node ids beneath the input node.
gather_leaves <- function(node, df) {
  leaves <- c()
  queue <- c(node)
  while (length(queue) > 0) {
    current <- queue[1]
    quelen <- length(queue)
    if (quelen == 1) {
      queue <- c()
    } else {
      queue <- queue[2:quelen]
    }
    
    children <- df$node[df$parent == current & !is.na(df$parent)]
    if (length(children) == 0) {
      leaves <- c(leaves, current)
    } else {
      queue <- c(queue, children)
    }
  }
  return(sort(leaves))
}

# Returns the direct children of a node that are themselves internal (have children).
# Input: node (scalar node id), df (data.frame with at least columns: node, parent).
# Output: sorted integer vector of "latent" (non-leaf) child node ids.
gather_direct_latent_nodes <- function(node, df) {
  latent_children <- c()
  children <- df$node[df$parent == node & !is.na(df$parent)]
  for (child in children) {
    grand_children <- df$node[df$parent == child & !is.na(df$parent)]
    if (length(grand_children) > 0) {
      latent_children <- c(latent_children, child)
    }
  }
  return(sort(latent_children))
}

# Builds per-depth lists of groups, weights, children, and statuses from a grouped node table.
# Input: group_df (data.frame with columns: node, depth, leaves(list), latent_children(list), weight).
# Output: list with elements: groups (list of integer vectors), weights (numeric lists), children (list), status (list).
make_list=function(group_df){
  weight_list=list()
  group_list=list()
  children_list=list()
  status_list=list()
  depth=sort(unique(group_df$depth),decreasing = T)
  for (i in 1:length(depth)) {
    temp_df=group_df[group_df$depth==depth[i],]
    
    status_list[[i]]=rep(1,length(temp_df$node))
    names(status_list[[i]])=temp_df$node
    
    group_list[[i]]=temp_df$leaves
    names(group_list[[i]])=temp_df$node
    
    weight_list[[i]]=temp_df$weight
    names(weight_list[[i]])=temp_df$node
    
    children_list[[i]]=(temp_df$latent_children)
    names(children_list[[i]])=temp_df$node
    
  }
  return(list(groups=group_list,weights=weight_list,children=children_list,status=status_list))
}

# For each non-leaf, gathers its leaf set, latent children, and weight; plus a packed list format.
# Input: df (data.frame with at least columns: node, parent, weight).
# Output: list(list_result = make_list(...) structure, df_result = data.frame with node, depth, leaves, latent_children, weight).
gather_leaf_nodes_per_non_leaf <- function(df) {
  non_leaf_nodes_df <- detect_non_leaf_nodes(df)
  non_leaf_nodes_df$leaves <- lapply(non_leaf_nodes_df$node, gather_leaves, df = df)
  non_leaf_nodes_df$latent_children <- lapply(non_leaf_nodes_df$node, gather_direct_latent_nodes, df = df)
  non_leaf_nodes_df$weight <- sapply(non_leaf_nodes_df$node, function(node) df$weight[df$node == node])
  result <- make_list(non_leaf_nodes_df)
  return(list(list_result = result, df_result = non_leaf_nodes_df))
}

############################################## This section includes algorithms for squares loss #################################################

# Accelerated proximal gradient (FISTA) solver for squared loss + ridge + tree-guided aggregation penalty (via prox_tree with backtracking).
# Input: Y (n-vector), X (n*p), tree_result (list from gather_leaf_nodes_per_non_leaf)$list_result, lambda (numeric), warm_start (bool), init_beta (p-vector), intercept (bool), ridge_param (numeric), thresh (numeric).
# Output: if intercept=FALSE, numeric p-vector of coefficients; else list(beta0 = intercept, beta1 = p-vector).
acc_prox_simple_linear=function(Y,X,tree_result,lambda,warm_start=F,init_beta=NULL,intercept=F,ridge_param=0,thresh=1e-6){
  if(intercept){
    Y.mean=mean(Y)
    X.mean=apply(X, 2, mean)
    Y=Y-Y.mean
    X=scale(X,scale = F)
  }
  
  p=ncol(X)
  n=nrow(X)
  if (warm_start & !is.null(init_beta)) {
    beta0 = beta1 = init_beta
  } else {
    beta0 = beta1 = rep(1, p)
  }
  
  alpha0=1
  alpha1=0.5
  #tao=1
  
  matXX=crossprod(X)
  L0=eigen(matXX)$values[1]/n
  tao=n/L0
  
  matXY=crossprod(X,Y)
  
  iter=0
  dis=1
  consecutive_below_thresh = 0
  
  while (consecutive_below_thresh < 5) {
    iter=iter+1
    accept=F
    Gam=beta0+(alpha0-1)/alpha1*(beta1-beta0)
    
    matXXGam=crossprod(matXX,Gam)
    matXGam=crossprod(t(X),Gam)
    base1=-matXY+matXXGam+ridge_param*Gam
    while (!accept){
      eta=Gam-tao*base1/n
      eta.new=prox_tree(eta,lambda = tao*lambda, tree_result)
      if(sum((Y- crossprod(t(X),eta.new))^2)/(2*n) <= sum((Y- matXGam)^2)/(2*n) + sum((matXXGam- matXY)*(eta.new-Gam)) + sum((eta.new-Gam)^2)/(2*tao) || tao <= 1/L0){
        accept=T
      }
      else{
        tao=max(tao/2,1/L0)
      }
    }
    beta0=beta1
    beta1=eta.new
    old.obj=sum((Y- crossprod(t(X),beta0))^2)/(2*n)
    new.obj=sum((Y- crossprod(t(X),beta1))^2)/(2*n)
    dis=abs((new.obj-old.obj)/old.obj)
    
    if (dis < thresh) {
      consecutive_below_thresh = consecutive_below_thresh + 1
    } else {
      consecutive_below_thresh = 0
    }
    
    alpha0=alpha1
    alpha1=(1+sqrt(1+4*alpha0^2))/2
  }
  
  
  if(!intercept){
    return(as.vector(beta1)) 
  }
  else{
    return(list(beta0=Y.mean-sum(X.mean*beta1),beta1=as.vector(beta1)))
  }
}

# Counts the total number of leaf nodes p in a (possibly multi-root) tree.
# Input: tree_df (data.frame with at least columns: node, parent).
# Output: integer p = total leaf count.
find_p=function(tree_df){
  roots=tree_df$node[is.na(tree_df$parent)]
  p=0
  for(i in 1:length(roots)){
    p=p+length(gather_leaves(roots[i],tree_df))
  }
  return(p)
}

# Builds the coarsest covering set of (weighted) groups from a tree to span all leaves, preferring shallow/nonzero-weight nodes.
# Input: df (original tree_df with node/parent/weight), result (non-leaf summary data.frame with columns: node, depth, leaves, latent_children, weight).
# Output: data.frame coarest_set with rows of selected nodes and columns (node, depth, leaves, latent_children, weight).
find_coarest=function(df, result){
  p=find_p(df)
  all_cover = F
  max_depth = max(result$depth)
  coarest_set = result[result$depth == 0, ]
  if (result[result$depth == 0, ]$weight != 0) {
    indiv.nodes=(1:p)[-sort(unlist((coarest_set$leaves)))]
    if(length(indiv.nodes)!=0){
      for(i in 1:length(indiv.nodes)){
        new_df=data.frame(node=indiv.nodes[i],depth=1, leaves=list(indiv.nodes[i]),latent_children=NA,weight=0)
        names(new_df)=names(coarest_set)
        #print(names(new_df))
        coarest_set=rbind(coarest_set,new_df)
      }
    }
    return(coarest_set)
  }
  current_depth = 1
  while (current_depth <= max_depth && all_cover == F) {
    sub_result = result[result$depth == current_depth, ]
    sub_result = sub_result[!is.element(df$parent[sub_result$node], 
                                        coarest_set$node[-1]), ]
    if (nrow(sub_result) == 0) {
      if (coarest_set$weight[1] == 0) 
        coarest_set = coarest_set[-1, ]
      
      indiv.nodes=(1:p)[-sort(unlist((coarest_set$leaves)))]
      if(length(indiv.nodes)!=0){
        for(i in 1:length(indiv.nodes)){
          new_df=data.frame(node=indiv.nodes[i],depth=1,leaves=list(indiv.nodes[i]),latent_children=NA,weight=0)
          
          names(new_df)=names(coarest_set)
          coarest_set=rbind(coarest_set,new_df)
        }
      }
      return(coarest_set)
    }
    coarest_set = rbind(coarest_set, sub_result[sub_result$weight != 
                                                  0, ])
    if (length(unique(unlist(coarest_set[-1, ]$leaves))) == 
        length(coarest_set[1, ]$leaves[[1]])) 
      all_cover = T
    current_depth = current_depth + 1
  }
  if (coarest_set$weight[1] == 0) 
    coarest_set = coarest_set[-1, ]
  
  indiv.nodes=(1:p)[-sort(unlist((coarest_set$leaves)))]
  if(length(indiv.nodes)!=0){
    for(i in 1:length(indiv.nodes)){
      new_df=data.frame(node=indiv.nodes[i],depth=1,leaves=list(indiv.nodes[i]),latent_children=NA,weight=0)
      names(new_df)=names(coarest_set)
      coarest_set=rbind(coarest_set,new_df)
    }
  }
  return(coarest_set)
}

# Computes lambda_max (penalty upper bound) for the tree penalty via grouped collapse Q, OLS on XQ, and KKT group scores.
# Input: Y (n-vector), X (n*p), coarest_set (data.frame from find_coarest with at least columns: leaves (list), weight).
# Output: numeric scalar penalty.max = max group score / n.
find_max_param.linear=function(Y,X,coarest_set){
  p=ncol(X)
  n=nrow(X)
  leaves=sort(coarest_set[1,]$leaves[[1]], decreasing = F)
  valid_set=coarest_set[coarest_set$weight!=0,]
  p1=p-sum(unlist(lapply(valid_set$leaves, length)))+nrow(valid_set)
  Q=matrix(0,nrow = p, ncol = p1 )
  remain=rep(1,length(leaves))
  for (i in 1:nrow(valid_set)) {
    ind=match(valid_set[i,]$leaves[[1]],leaves)
    Q[ind,i]=1
    remain[ind]=0
  }
  #print(remain)
  indiv_num=sum(remain)
  if(indiv_num>0) Q[which(remain==1),(p1- indiv_num+1):p1]=diag(indiv_num)
  #print(Q)
  X1=X%*%Q
  #print(X1)
  beta1=ginv(t(X1)%*%X1)%*%t(X1)%*%Y
  target=t(X)%*%(Y-X1%*%beta1)
  #print(target)
  vals=sqrt((t(Q)%*%target^2)[1:nrow(valid_set)])/valid_set$weight
  return(max(vals)/n)
}

# Solves along a lambda grid (warm starts) for the tree-guided model and records path and coefficient MSE vs true_beta.
# Input: Y (n-vector), X (n*p), tree_df (with node/parent/weight), true_beta (p-vector), ridge.param (numeric), seqc (optional lambda vector), thresh (numeric).
# Output: list(loss = length(seqc) vector of MSE vs true_beta, beta = p*length(seqc) matrix of coefficients, lambda = lambda grid).
grid.simple_linear=function(Y,X,tree_df,true_beta,ridge.param=0,seqc=NULL,thresh=1e-5){
  n=length(Y)
  p=ncol(X)
  tree_result=gather_leaf_nodes_per_non_leaf(tree_df)
  coarest_set=find_coarest(tree_df,tree_result$df_result)
  penalty.max=find_max_param.linear(Y-mean(Y),scale(X,scale = F),coarest_set)
  if(is.null(seqc)){
    seqc=exp(seq(-4-log(n),log(penalty.max/5),length=50))
  }
  Y1=Y
  X1=X
  beta=matrix(0,nrow = p ,ncol = length(seqc))
  beta[,1]=acc_prox_simple_linear(Y1, X1,tree_result$list_result,lambda = seqc[1],ridge_param = ridge.param,thresh = thresh)
  for (i in 2:length(seqc)) {
    beta[,i]=acc_prox_simple_linear(Y1, X1,tree_result$list_result,lambda = seqc[i],warm_start = T,init_beta = as.vector(beta[,i-1]),ridge_param = ridge.param,thresh = thresh)
  }
  loss=apply(matrix(rep(true_beta,length(seqc)),nrow = p,byrow = F)-beta,2,function(x) sum(x^2))/p
  return(list(loss=loss,beta=beta,lambda=seqc))
}


# Cross-validates the tree-guided squared-loss model to pick lambda, then refits on all data.
# Input: Y (n), X (n*p), tree_df (tree metadata), folds, seqc, thresh, ridge.param, intercept, Mmratio.
# Output: list(selected.param, beta, valid.error).
cv.simple_linear=function (Y, X, tree_df, folds = 5, seqc = NULL,
    thresh = 1e-05, ridge.param = 0, intercept = F, Mmratio = 10000) 
{
    n = length(Y)
    p = ncol(X)
    stopifnot(n >= 2 * folds)
    
    random_sequence <- sample(1:n)
    index <- cut(random_sequence, breaks = folds, labels = FALSE)
    tree_result = gather_leaf_nodes_per_non_leaf(tree_df)
    coarest_set = find_coarest(tree_df, tree_result$df_result)
    penalty.max = find_max_param.linear(Y, X, coarest_set)
    if (is.null(seqc)) {
        seqc = exp(seq(log(penalty.max/Mmratio), log(penalty.max), 
            length = 50))
    }
    seqc = sort(seqc, decreasing = T)
    vals.mat = matrix(0, nrow = folds, ncol = length(seqc))
    colnames(vals.mat) = seqc
    for (i in 1:folds) {
        print(i)
        X_train = X[which(index != i), ]
        X_test = X[which(index == i), ]
        Y_train = Y[which(index != i)]
        Y_test = Y[which(index == i)]
        res = grid.simple_linear(Y = Y_train, X = X_train, tree_df = tree_df, 
            true_beta = rep(0, p), seqc = seqc, thresh = thresh, 
            ridge.param = ridge.param)
        #vals.mat[i, ] = negtv_lglkh(Y_test, X_test, beta = res$beta)
        vals.mat[i, ] = apply(crossprod(t(X_test),res$beta),2,function(x) mean((Y_test-x)^2))
    }
    vals.vec = apply(vals.mat, 2, mean)
    selected.param = seqc[which.min(vals.vec)]
    beta = acc_prox_simple_linear(Y, X, tree_result$list_result, 
        lambda = selected.param, intercept = intercept, ridge_param = ridge.param)
    return(list(selected.param = selected.param, beta = beta, 
        valid.error = vals.vec))
}

                              
# Cross-validates RARE (linear/Gaussian) to pick lambda, then refits on all data.
# Input: Y (n), X (n*p), tree_df, folds, seqc, thresh, ridge.param, intercept.
# Output: list(selected.param, beta, valid.error).
cv.rare.linear=function (Y, X, tree_df, folds = 5, seqc = NULL, thresh = 1e-05, ridge.param = 0, intercept = F) 
{
    n = length(Y)
    p = ncol(X)
    stopifnot(n >= 2 * folds)
    random_sequence <- sample(1:n)
    index <- cut(random_sequence, breaks = folds, labels = FALSE)
    A=df_to_A(tree_df = tree_df,p)
    A_sparse <- as(A, "dgCMatrix")
    if (is.null(seqc)) {
        seqc = rarefit(y = Y, X = X, A=A_sparse, 
            intercept = F,alpha=1)$lambda
    }
    seqc = sort(seqc, decreasing = T)
    vals.mat = matrix(0, nrow = folds, ncol = length(seqc))
    colnames(vals.mat) = seqc
    for (i in 1:folds) {
        print(i)
        X_train = X[which(index != i), ]
        X_test = X[which(index == i), ]
        Y_train = Y[which(index != i)]
        Y_test = Y[which(index == i)]
        res = rarefit(y = Y_train, X = X_train, A_sparse, intercept = F, lambda = seqc,alpha = 1)
        #print(length(res$beta))
        vals.mat[i, ] = apply(res$beta[[1]],2,function(x) mean((Y_test-crossprod(t(X_test),x))^2))
    }
    vals.vec = apply(vals.mat, 2, mean)
    selected.param = seqc[which.min(vals.vec)]
    beta = as.numeric(rarefit(Y, X, A_sparse, intercept = F, lambda = selected.param,alpha=1)$beta[[1]])
    return(list(selected.param = selected.param, beta = beta, 
        valid.error = vals.vec))
}


# Estimates an in-sample signal-to-noise ratio for RARE via K-fold CV.
# Input: y (n), X (n*p), tree_df, folds.
# Output: mean SNR across folds (numeric scalar).
INratio_rare=function(y,X,tree_df,folds=5){
    n=nrow(X)
    p=ncol(X)
    random_sequence <- sample(1:n)
    index <- cut(random_sequence, breaks = folds, labels = FALSE)
    snr_seq=numeric(folds)
    for(i in 1:folds){
        X.train=X[index!=i,]
        y.train=y[index!=i]
        X.test=X[index==i,]
        y.test=y[index==i]
        #print(dim(X.train))
        #print(length(y.train))
        
        model=cv.rare.linear(y.train,X.train,tree_df,folds = folds,intercept = F)
        beta=model$beta
        #print(beta)
        predictions <- as.vector(X.test %*% beta)
        mse <- mean((y.test - predictions)^2)
        var_pred <- var(predictions)
        snr_seq[i] <- var_pred / mse
    }
    return(mean(snr_seq))
}

                              
# One train/validation split comparing LASSO, Ridge, RARE, and our tree-guided model (Gaussian loss).
# Input: y (n), X (n*p), tree_df, split.seed, nfolds, nval, ridge_param, Mmratio.
# Output: list of coefficients, losses, and selected lambda for each method (+ optional “new LASSO” on grouped X).
real_data_one_round.linear=function (y, X, tree_df, split.seed = 123, nfolds = 5, nval = 31, 
    ridge_param = 0, Mmratio = 10000) 
{
    n = nrow(X)
    p = ncol(X)
    set.seed(split.seed)
    index = sample(1:nrow(X), nval)
    X.train = X[-index, ]
    X.test = X[index, ]
    y.train = y[-index]
    y.test = y[index]
    set.seed(2 * split.seed)
    lasso.mod = cv.glmnet(X.train, y.train, family = "gaussian", 
        alpha = 1, intercept = F, nfolds = nfolds)
    lasso.beta = as.numeric(glmnet(X.train, y.train, family = "gaussian", 
        alpha = 1, intercept = F, lambda = lasso.mod$lambda.min)$beta)
    lasso.loss = mean((y.test - crossprod(t(X.test), lasso.beta))^2)
    lasso.param = lasso.mod$lambda.min
    set.seed(2 * split.seed)
    ridge.mod = cv.glmnet(X.train, y.train, family = "gaussian", 
        alpha = 0, intercept = F, nfolds = nfolds)
    ridge.beta = as.numeric(glmnet(X.train, y.train, family = "gaussian", 
        alpha = 0, intercept = F, lambda = ridge.mod$lambda.min)$beta)
    ridge.loss = mean((y.test - crossprod(t(X.test), ridge.beta))^2)
    ridge.param = ridge.mod$lambda.min
    set.seed(2 * split.seed)
    rare.mod = cv.rare.linear(y.train, X.train, tree_df, folds = nfolds)
    rare.beta = as.numeric(rare.mod$beta)
    rare.loss = mean((y.test - crossprod(t(X.test), rare.beta))^2)
    rare.param = rare.mod$selected.param
    set.seed(2 * split.seed)
    our.mod = cv.simple_linear(y.train, X.train, tree_df, folds = nfolds, 
        ridge.param = ridge_param, Mmratio = Mmratio)
    our.beta = our.mod$beta
    our.loss = mean((y.test - crossprod(t(X.test), our.beta))^2)
    our.param = our.mod$selected.param
    
    if(length(table(our.beta))>1){
        our.levels=factor(our.beta,labels = 1:length(table(our.beta)))
        new.X.train=matrix(0,nrow = nrow(X.train),ncol=length(table(our.beta)))
        new.X.test=matrix(0,nrow = nrow(X.test),ncol=length(table(our.beta)))
        for(i in 1:length(table(our.beta))){
            new.X.train[,i]=rowSums(as.matrix(X.train[,which(our.levels==i)],nrow=nrow(X.train)))
            new.X.test[,i]=rowSums(as.matrix(X.test[,which(our.levels==i)],nrow=nrow(X.test)))
        }
        set.seed(2 * split.seed)
        new.lasso.mod=cv.glmnet(new.X.train, y.train, family = "gaussian", 
        alpha = 1, intercept = F, nfolds = nfolds)
        new.lasso.beta = as.numeric(glmnet(new.X.train, y.train, family = "gaussian", 
        alpha = 1, intercept = F, lambda = lasso.mod$lambda.min)$beta)
        new.lasso.loss = mean((y.test - crossprod(t(new.X.test), new.lasso.beta))^2)
        new.lasso.param = new.lasso.mod$lambda.min
    }else{
        new.lasso.beta=new.lasso.loss=new.lasso.param=NA
    }
    return(list(lasso.beta = lasso.beta, lasso.loss = lasso.loss, 
        lasso.param = lasso.param, ridge.beta = ridge.beta, ridge.loss = ridge.loss, 
        ridge.param = ridge.param, rare.beta = rare.beta, rare.loss = rare.loss, 
        rare.param = rare.param, our.beta = our.beta, our.loss = our.loss, 
        our.param = our.param,new.lasso.beta = new.lasso.beta, new.lasso.loss = new.lasso.loss, 
        new.lasso.param = new.lasso.param))
}


###############################################  This section includes algorithms for entropy loss ################################################

#FISTA solver for logistic loss with ridge and tree-guided proximal step (with backtracking and optional warm start).
#Input: Y (n), X (n*p), tree_result (from gather_leaf_nodes_per_non_leaf)$list_result, lambda, warm_start, init_beta, intercept, ridge_param, thresh.
#Output: beta (p) if intercept=FALSE; otherwise list(beta0, beta1).
acc_prox_simple_logistic=function(Y,X,tree_result,lambda,warm_start=F,init_beta=NULL,intercept=F,ridge_param=0,thresh=1e-6){
  if(intercept){
    Y.mean=mean(Y)
    X.mean=apply(X, 2, mean)
    Y=Y-Y.mean
    X=scale(X,scale = F)
  }
  
  p=ncol(X)
  n=nrow(X)
  if (warm_start & !is.null(init_beta)) {
    beta0 = beta1 = init_beta
  } else {
    beta0 = beta1 = rep(1, p)
  }
  
  alpha0=1
  alpha1=0.5
  
  matXX=crossprod(X)
  L0=eigen(matXX)$values[1]/n
  tao=n/L0
  
  matXY=crossprod(X,Y)
  
  iter=0
  dis=1
  consecutive_below_thresh = 0
  
  while (consecutive_below_thresh < 5 && iter<1e6) {
    iter=iter+1
    accept=F
    Gam=beta0+(alpha0-1)/alpha1*(beta1-beta0)
    
    matXGam=crossprod(t(X),Gam)
    exp.matXGam=exp(matXGam)
    base1=crossprod(X,exp.matXGam/(1+exp.matXGam)-Y)+ridge_param*Gam
    while (!accept){
      eta=Gam-tao*base1/n
      eta.new=prox_tree(eta,lambda = tao*lambda, tree_result)
      if(mean(log(1+exp(crossprod(t(X),eta.new))))-mean(crossprod(t(X),eta.new)*Y) <= mean(log(1+exp.matXGam))-mean(matXGam*Y) + sum((crossprod(X,exp.matXGam/(1+exp.matXGam)-Y)/n)*(eta.new-Gam)) + sum((eta.new-Gam)^2)/(2*tao) || tao <= 1/L0){
        accept=T
      }
      else{
        tao=max(tao/2,1/L0)
      }
    }
    beta0=beta1
    beta1=eta.new
    old.obj=mean(log(1+exp(crossprod(t(X),beta0))))-mean(crossprod(t(X),beta0)*Y)
    new.obj=mean(log(1+exp(crossprod(t(X),beta1))))-mean(crossprod(t(X),beta1)*Y)
    dis=abs((new.obj-old.obj)/old.obj)
    
    if (dis < thresh) {
      consecutive_below_thresh = consecutive_below_thresh + 1
    } else {
      consecutive_below_thresh = 0
    }
    
    alpha0=alpha1
    alpha1=(1+sqrt(1+4*alpha0^2))/2
  }
  
  if(!intercept){
    return(as.vector(beta1)) 
  }
  else{
    return(list(beta0=Y.mean-sum(X.mean*beta1),beta1=as.vector(beta1)))
  }
}

#Computes lambda_max for the tree penalty via grouped collapse Q and a logistic GLM fit on XQ.
#Input: Y (n), X (n*p), coarest_set (from find_coarest with leaves/weight).
#Output: numeric scalar lambda_max/n.
find_max_param.logistic=function (Y, X, coarest_set) {
    p = ncol(X)
    n = nrow(X)
    leaves = 1:p
    valid_set = coarest_set[coarest_set$weight != 0, ]
    p1 = p - sum(unlist(lapply(valid_set$leaves, length))) + 
        nrow(valid_set)
    Q = matrix(0, nrow = p, ncol = p1)
    remain = rep(1, length(leaves))
    for (i in 1:nrow(valid_set)) {
        ind = match(valid_set[i, ]$leaves[[1]], leaves)
        Q[ind, i] = 1
        remain[ind] = 0
    }
    indiv_num = sum(remain)
    if (indiv_num > 0) 
        Q[which(remain == 1), (p1 - indiv_num + 1):p1] = diag(indiv_num)
    X1 = X %*% Q
    fit <- glm(Y ~ X1 + 0, family = binomial(link = "logit"))
    beta1 = coef(fit)
    target = crossprod(X, Y - 1/(1 + exp(-crossprod(t(X1), beta1))))
    vals = sqrt((t(Q) %*% target^2)[1:nrow(valid_set)])/valid_set$weight
    return(max(vals)/n)
}


#Solves the tree-guided logistic model along a lambda-grid (warm starts) and records coefficient paths and MSE vs true beta.
#Input: Y (n), X (n*p), tree_df, true_beta (p), ridge.param, seqc (optional lambda vector), thresh.
#Output: list(loss (|lambda|), beta (p*|lambda|), lambda (|lambda|)).
grid.logistic=function(Y,X,tree_df,true_beta,ridge.param=0,seqc=NULL,thresh=1e-5){
  n=length(Y)
  p=ncol(X)
  tree_result=gather_leaf_nodes_per_non_leaf(tree_df)
  coarest_set=find_coarest(tree_df,tree_result$df_result)
  penalty.max=find_max_param.logistic(Y,X,coarest_set)
  if(is.null(seqc)){
    seqc=exp(seq(-4-log(n),log(penalty.max),length=50))
  }
  Y1=Y
  X1=X
  beta=matrix(0,nrow = p ,ncol = length(seqc))
  beta[,1]=acc_prox_simple_logistic(Y1, X1,tree_result$list_result,lambda = seqc[1],ridge_param = ridge.param,thresh = thresh)
  for (i in 2:length(seqc)) {
    beta[,i]=acc_prox_simple_logistic(Y1, X1,tree_result$list_result,lambda = seqc[i],warm_start = T,init_beta = as.vector(beta[,i-1]),ridge_param = ridge.param,thresh = thresh)
  }
  loss=apply(matrix(rep(true_beta,length(seqc)),nrow = p,byrow = F)-beta,2,function(x) sum(x^2))/p
  return(list(loss=loss,beta=beta,lambda=seqc))
}
             

#Fits logistic RARE with a tree-expanded design A (built from tree_df) and group-wise penalty factors.
#Input: y (n), X (n*p), tree_df, A (ignored; built internally), hc (unused), intercept, lambda (optional), nlam, lam.min.ratio, eps, maxite.
#Output: list(beta0 (0s), beta (p×|lambda|), gamma (|A|×|lambda|), lambda, A, intercept).
rarefit.logistic=function (y, X, tree_df = NULL,A=NULL,hc, intercept = F, lambda = NULL, nlam = 50, 
    lam.min.ratio = 1e-04, eps = 1e-05, maxite = 1e+06) 
{
    n <- nrow(X)
    p <- ncol(X)
    A=df_to_A(tree_df,p)
    tree_result=gather_leaf_nodes_per_non_leaf(tree_df)$df_result
    tree_result=tree_result[order(tree_result$depth,decreasing = T), ]
    zero.ind=which(tree_result$depth==0)+p
    nnodes <- ncol(A)
    penalty.factor=rep(1,nnodes)
    penalty.factor[zero.ind]=0
    X_use <- X <- as.matrix(X)
    y_use <- as.vector(y)
    if (is.null(lambda)) {
        lambda <- max(abs(t(X_use %*% A) %*% (y_use - 1/2)))/n * 
            exp(seq(0, log(lam.min.ratio), len = nlam))
    }
    else {
        if (min(lambda) < 0) 
            stop("lambda cannot be negative.")
        nlam <- length(lambda)
    }
    beta0 = numeric(nlam)
    beta <- gamma <- c()
    ret <- glmnet(X_use %*% A, y_use, family = "binomial", lambda = lambda, 
        standardize = F, intercept = F, penalty.factor = penalty.factor, thresh = eps, maxit = maxite)
    beta <- as.matrix(A %*% ret$beta)
    gamma <- as.matrix(ret$beta)
    list(beta0 = beta0, beta = beta, gamma = gamma, lambda = lambda, 
        A = A, intercept = intercept)
}

             

#Computes average negative log-likelihood for logistic regression at given coefficients.
#Input: Y (n or scalar), X (n*p or 1*p), beta (p*m or p).
#Output: numeric (or length-m numeric) negative log-likelihood.
negtv_lglkh=function(Y,X,beta){
    n=nrow(X)
    temp=crossprod(t(X),beta)
    return(1/n*as.numeric(-crossprod(Y,temp)+colSums(log(1+exp(temp)))))
}
             

#K-fold CV to select lambda for the tree-guided logistic model; returns selected lambda and refit coefficients.
#Input: Y (n), X (n*p), tree_df, folds, seqc (optional lambda grid), neg.ind (optional indices for stratification), thresh, ridge.param, intercept, Mmratio.
#Output: list(selected.param, beta (p), valid.error (per-lambda mean loss)).
cv.logistic=function(Y,X,tree_df,folds=5,seqc=NULL,neg.ind=NULL,thresh=1e-5,ridge.param=0,intercept=F,Mmratio=1e+4){
  n=length(Y)
  p=ncol(X)
  stopifnot(n>=2*folds)

  if(is.null(neg.ind)){
      random_sequence <- sample(1:n)
      index <- cut(random_sequence, breaks=folds, labels=FALSE)
  }else{
      pos.ind.rand=sample((1:n)[-neg.ind])
      neg.ind.rand=sample(neg.ind)
      pos.partition=cut(pos.ind.rand, breaks=folds, labels=FALSE)
      neg.partition=cut(neg.ind.rand, breaks=folds, labels=FALSE)
      index=numeric(n)
      index[-neg.ind]=pos.partition
      index[neg.ind]=neg.partition
  }
  #print(index)
    
    
  tree_result=gather_leaf_nodes_per_non_leaf(tree_df)
  coarest_set=find_coarest(tree_df,tree_result$df_result)
  penalty.max=find_max_param.logistic(Y,X,coarest_set)
  #print(penalty.max)
  if(is.null(seqc)){
    seqc=exp(seq(log(penalty.max/Mmratio),log(penalty.max),length=50))
  }
  seqc=sort(seqc,decreasing = T)
  vals.mat=matrix(0,nrow = folds,ncol = length(seqc))
  colnames(vals.mat)=seqc
  
  for(i in 1:folds){
    print(i)
    X_train=X[which(index!=i),]
    X_test=X[which(index==i),]
    
    Y_train=Y[which(index!=i)]
    Y_test=Y[which(index==i)]
      
    res=grid.logistic(Y = Y_train,X = X_train,tree_df = tree_df,true_beta = rep(0,p),seqc = seqc,thresh = thresh,ridge.param = ridge.param)
    vals.mat[i,]=negtv_lglkh(Y_test,X_test,beta = res$beta)
  }
    #print(vals.mat)
  vals.vec=apply(vals.mat, 2, mean)
  selected.param=seqc[which.min(vals.vec)]
  beta=acc_prox_simple_logistic(Y,X,tree_result$list_result,lambda = selected.param,intercept = intercept,ridge_param = ridge.param)
  return(list(selected.param=selected.param,beta=beta,valid.error=vals.vec))
}

#Leave-one-out CV over a lambda-grid for the tree-guided logistic model.
#Input: Y (n), X (n*p), tree_df, seqc (optional lambda grid).
#Output: list(cve (n*|lambda| losses), selected.param (length n), seqc).
loocv.logistic=function(Y,X,tree_df,seqc=NULL){
    n=length(Y)
    p=ncol(X)
    
    tree_result=gather_leaf_nodes_per_non_leaf(tree_df)
  coarest_set=find_coarest(tree_df,tree_result$df_result)
  penalty.max=find_max_param.logistic(Y,X,coarest_set)
  #print(penalty.max)
  if(is.null(seqc)){
    seqc=exp(seq(log(penalty.max/1e+4),log(penalty.max),length=50))
  }
    cve=matrix(0,nrow=n,ncol=length(seqc))
    selected.param=numeric(n)
    for(i in 1:n){
        print(i)
        res=grid.logistic(Y[-i],X[-i,],tree_df,true_beta = rep(0,p),seqc = seqc,ridge.param = ridge.param)
        loss=negtv_lglkh(Y[i],matrix(X[i,],nrow=1),res$beta)
        #loss=(1/(1+exp(-crossprod(X[i,],res$beta)))>0.5)!=Y[i]
        cve[i,]=loss
        selected.param[i]=which.min(loss)
    }
    
    return(list(cve=cve,selected.param=selected.param,seqc=seqc))
}


#Builds an expanded design mapping matrix A from tree_df (columns: p leaves + one column per nonzero-weight internal node).
#Input: tree_df, p (number of leaves/features).
#Output: A (p*(p+g)) 0/1 matrix mapping features to groups.
df_to_A=function (tree_df, p){
    tree_result = gather_leaf_nodes_per_non_leaf(tree_df)$df_result
    tree_result=tree_result[order(tree_result$depth,decreasing = T), ]
    A = diag(p)
    for (i in 1:nrow(tree_result)) {
        new.col = rep(0, p)
        if (tree_result[i, "weight"] != 0) {
            new.col[tree_result[i, "leaves"][[1]]] = 1
            A = cbind(A, new.col)
        }
    }
    return(A)
}

#K-fold CV for RARE (logistic) over a supplied lambda grid (or one inferred from data).
#Input: Y (n), X (n*p), tree_df, folds, seqc (optional), neg.ind, thresh, ridge.param, intercept.
#Output: list(selected.param, beta (rarefit object), valid.error).
cv.rare.logistic=function (Y, X, tree_df, folds = 5, seqc = NULL, neg.ind = NULL, 
    thresh = 1e-05, ridge.param = 0, intercept = F) 
{
    n = length(Y)
    p = ncol(X)
    stopifnot(n >= 2 * folds)
    if (is.null(neg.ind)) {
        random_sequence <- sample(1:n)
        index <- cut(random_sequence, breaks = folds, labels = FALSE)
    }
    else {
        pos.ind.rand = sample((1:n)[-neg.ind])
        neg.ind.rand = sample(neg.ind)
        pos.partition = cut(pos.ind.rand, breaks = folds, labels = FALSE)
        neg.partition = cut(neg.ind.rand, breaks = folds, labels = FALSE)
        index = numeric(n)
        index[-neg.ind] = pos.partition
        index[neg.ind] = neg.partition
    }
    if (is.null(seqc)) {
        seqc = rarefit.logistic(y = Y, X = X, tree_df=tree_df, intercept = F)$lambda
    }
    seqc = sort(seqc, decreasing = T)
    vals.mat = matrix(0, nrow = folds, ncol = length(seqc))
    colnames(vals.mat) = seqc
    for (i in 1:folds) {
        print(i)
        X_train = X[which(index != i), ]
        X_test = X[which(index == i), ]
        Y_train = Y[which(index != i)]
        Y_test = Y[which(index == i)]
        res = rarefit.logistic(y = Y_train, X = X_train, tree_df=tree_df, 
            intercept = F, lambda = seqc)
        vals.mat[i, ] = negtv_lglkh(Y_test, X_test, beta = res$beta)
    }
    vals.vec = apply(vals.mat, 2, mean)
    selected.param = seqc[which.min(vals.vec)]
    beta = rarefit.logistic(Y, X, tree_df=tree_df, intercept = F, lambda = selected.param)
    return(list(selected.param = selected.param, beta = beta, 
        valid.error = vals.vec))
}

#Leave-one-out CV for RARE (logistic) across a lambda grid.
#Input: Y (n), X (n*p), tree_df, seqc (optional lambda grid).
#Output: list(cve (n*|lambda|), selected.param (length n), seqc).
loocv.rare.logistic=function(Y,X,tree_df,seqc=NULL){
    n=length(Y)
    p=ncol(X)
    A=df_to_A(tree_df,p)
    
    if(is.null(seqc)){
    seqc=rarefit.logistic(y = Y,X=X,A = A,intercept = F)$lambda
  }
  seqc=sort(seqc,decreasing = T)
    
    cve=matrix(0,nrow=n,ncol=length(seqc))
    selected.param=numeric(n)
    for(i in 1:n){
        res=rarefit.logistic(y=Y[-i],X=X[-i,],A = A,intercept = F,lambda = seqc)
        loss=negtv_lglkh(Y[i],matrix(X[i,],nrow=1),res$beta)
        #loss=(1/(1+exp(-crossprod(X[i,],res$beta)))>0.5)!=Y[i]
        cve[i,]=loss
        selected.param[i]=which.min(loss)
    }
    
    return(list(cve=cve,selected.param=selected.param,seqc=seqc))
}


#One split evaluation reporting loss and AUC for LASSO, Ridge, RARE, and the tree-guided logistic model.
#Input: y (n), X (n*p), tree_df, split.seed, nfolds, nval (test size), ridge_param, Mmratio.
#Output: list of each method’s beta, loss, AUC, and selected tuning parameter.
real_data_one_round_auc=function (y, X, tree_df, split.seed = 123, nfolds = 5,nval=31,ridge_param=0,Mmratio=1e+4) 
{
    n = nrow(X)
    p = ncol(X)
    set.seed(split.seed)
    index = sample(1:nrow(X), nval)
    X.train = X[-index, ]
    X.test = X[index, ]
    y.train = y[-index]
    y.test = y[index]
    set.seed(2 * split.seed)
    lasso.mod = cv.glmnet(X.train, y.train, family = "binomial", 
        alpha = 1, intercept = F, nfolds = nfolds)
    lasso.beta = as.numeric(glmnet(X.train, y.train, family = "binomial", 
        alpha = 1, intercept = F, lambda = lasso.mod$lambda.min)$beta)
    lasso.loss = negtv_lglkh(y.test, X.test, lasso.beta)
    lasso.auc = auc(roc(predictor = 1/(1 + exp(-crossprod(t(X.test), lasso.beta))),response = y.test)) 
    lasso.param = lasso.mod$lambda.min
    set.seed(2 * split.seed)
    ridge.mod = cv.glmnet(X.train, y.train, family = "binomial", 
        alpha = 0, intercept = F, nfolds = nfolds)
    ridge.beta = as.numeric(glmnet(X.train, y.train, family = "binomial", 
        alpha = 0, intercept = F, lambda = ridge.mod$lambda.min)$beta)
    ridge.loss = negtv_lglkh(y.test, X.test, ridge.beta)
    ridge.auc = auc(roc(predictor = 1/(1 + exp(-crossprod(t(X.test), ridge.beta))),response = y.test)) 
    ridge.param = ridge.mod$lambda.min
    set.seed(2 * split.seed)
    rare.mod = cv.rare.logistic(y.train, X.train, tree_df, folds = nfolds)
    rare.beta = as.numeric(rare.mod$beta$beta)
    rare.loss = negtv_lglkh(y.test, X.test, rare.beta)
    rare.auc = auc(roc(predictor = 1/(1 + exp(-crossprod(t(X.test), rare.beta))),response = y.test)) 
    rare.param = rare.mod$selected.param
    set.seed(2 * split.seed)
    our.mod = cv.logistic(y.train, X.train, tree_df, folds = nfolds,ridge.param=ridge_param,Mmratio=Mmratio)
    our.beta = our.mod$beta
    our.loss = negtv_lglkh(y.test, X.test, our.beta)
    our.auc = auc(roc(predictor = 1/(1 + exp(-crossprod(t(X.test), our.beta))),response = y.test)) 
    our.param = our.mod$selected.param
    return(list(lasso.beta = lasso.beta, lasso.loss = lasso.loss, 
        lasso.auc = lasso.auc, lasso.param = lasso.param, ridge.beta = ridge.beta, 
        ridge.loss = ridge.loss, ridge.auc = ridge.auc, ridge.param = ridge.param, 
        rare.beta = rare.beta, rare.loss = rare.loss, rare.auc = rare.auc, 
        rare.param = rare.param, our.beta = our.beta, our.loss = our.loss, 
        our.auc = our.auc, our.param = our.param))
}


################################################## This section simulates the data with continuous outcome ################################################################
# Simulates a k-group structure over p leaves by drawing latent locations and building an hclust tree.
# Input: k (number of groups), p (number of leaves), tao (sd scale for within-group latent noise).
# Output: list(tree = hclust object, group = integer vector of group labels for each leaf).
simulate_tree=function(k,p,tao=0.05){
  means=1/(1:k)
  group.index=c(cut(1:ceiling(p*3/4),breaks=k/2,labels=F),cut((ceiling(p*3/4)+1):p,breaks=k/2,labels=F)+k/2)
  group.size=table((group.index))
  latent=c()
  for (i in 1:k) {
    near.index=ifelse(i<k,i+1,i-1)
    latent=c(latent,rnorm(group.size[i],means[i],tao*abs(means[i]-means[near.index])))
  }
  tree=hclust(dist(latent))
  return(list(tree=tree,group=group.index))
}

# Converts an hclust object to a tidy node/parent data frame and assigns node weights.
# Input: hc (hclust), weight.order (numeric exponent for size-based weighting; e.g., -1 downweights large groups).
# Output: data.frame with columns {node, parent, name, weight}, one row per leaf/latent node.
hclust_to_df <- function(hc,weight.order=-1) {
  tree_df <- data.frame(
    node = integer(),
    parent = integer(),
    name = character(),
    weight = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Function to recursively traverse the tree and fill the data frame
  traverse <- function(node, parent, name_prefix, depth) {
    if (node < 0) {
      # If node is a leaf, use its original index
      node_id <- -node
      name <- paste("leaf", node_id, sep = "")
      weight <- 0  # You can adjust this as needed
    } else {
      # If node is a latent node, assign a new ID
      node_id <- nrow(hc$merge) + node+1
      name <- paste("latent", node_id, sep = "")
      weight <- depth * 10  # Adjust this as needed
    }
    
    # Add current node to data frame
    tree_df <<- rbind(tree_df, data.frame(
      node = node_id,
      parent = parent,
      name = name,
      weight = weight,
      stringsAsFactors = FALSE
    ))
    
    # Recursively traverse left and right children if they exist
    if (node > 0) {
      left_child <- hc$merge[node, 1]
      right_child <- hc$merge[node, 2]
      traverse(left_child, node_id, name, depth + 1)
      traverse(right_child, node_id, name, depth + 1)
    }
  }
  
  # Start traversal from the root
  traverse(nrow(hc$merge), NA, "root", 0)
  
  # Reorder data frame by node ID
  tree_df <- tree_df %>% arrange(node)
  tree_df=assign_weight(tree_df,length(hc$order) ,weight.order=weight.order)
  return(tree_df)
}

# Assigns per-node weights based on group sizes extracted from the tree structure.
# Input: tree_df (data.frame with node/parent/weight), p (number of leaves), weight.order (exponent for |leaves|^weight.order).
# Output: data.frame tree_df with updated, mean-normalized weights for non-leaf nodes.
assign_weight=function(tree_df,p,weight.order=-1){
  tree_df$weight[(p+1):nrow(tree_df)]=1
  nodes=tree_df$node[tree_df$weight!=0]
  result=gather_leaf_nodes_per_non_leaf(tree_df)$df_result
  depth=c(nrow(result))
  for (i in 1:nrow(result)) {
    depth[i]=(length(result$leaves[[i]]))^weight.order
    tree_df$weight[tree_df$node==result$node[i]]=depth[i]
  }
  tree_df$weight=tree_df$weight/mean(depth)
  return(tree_df)
}

# Simulates grouped linear-model data with optional X scaling and two SNR controls (global-energy or variance-based).
# Input: n (#samples), group.index (length-p labels), s (group sparsity), beta.pre (optional k-vector), ratio (SNR via sum(means^2)), new.ratio (SNR via var(means)), scale.X (bool).
# Output: list(X = n*p matrix, Y = n-vector, true.beta = p-vector (adjusted if scale.X), A = p*k group-membership matrix).
simulate_data=function(n,group.index,s=0,beta.pre=NULL,ratio=5,fix_sigma=NULL,new.ratio=NULL,scale.X=FALSE){
  k=length(unique(group.index))
  p=length(group.index)
  A=matrix(0,nrow = p,k)
  for (i in 1:k) {
    A[,i]=(group.index==i)
  }
  if(is.null(beta.pre)){
    beta0=runif(k,1.5,2.5)
  }else{
    beta0=beta.pre
  }
  nonzero.ind=sample(1:k,ceiling(k*(1-s)))
  nonzero=rep(0,k)
  nonzero[nonzero.ind]=rep(c(-1,1),k%/%2)[1:length(nonzero.ind)]
  beta0=beta0*nonzero
  beta=A%*%beta0
  #X=matrix(rpois(n*p,0.02 * 5 / (p/k)),nrow = n,ncol = p)
  X=matrix(rpois(n*p,0.02),nrow = n,ncol = p)
  if(scale.X){
    scale.factor=sqrt(eigen(crossprod(X))$values[1])
    X=X/scale.factor
    beta=beta*scale.factor
  }
  means=X%*%beta
  if(is.null(new.ratio)){
    sigma=sqrt(sum(means^2)/ratio/n)
  }else{
    sigma=sqrt(var(means)/new.ratio)
  }
  #print(sigma)
  if(!is.null(fix_sigma)) sigma=fix_sigma
  Y=means+rnorm(n,0,sigma)
  #X=X/sqrt(eigen(crossprod(X))$values[1])
  return(list(X=X,Y=Y,true.beta=beta,A=A))
}

# Builds a grouped design by aggregating columns of X according to tree_group and returns the mapping matrix.
# Input: X (n×p matrix), tree_group (length-p vector of group labels).
# Output: list(X1 = n*G grouped design with G=number of unique groups, Q = p*G 0/1 membership matrix).
make_new_X=function(X,tree_group){
  index=tree_group
  levels=unique(index)
  Q=matrix(0,nrow = length(tree_group),ncol = length(levels))
  for (i in 1:length(levels)) {
    Q[index==levels[i],i]=1
  }
  return(list(X1=X%*%Q,Q=Q))
}

# Computes a Moore–Penrose pseudoinverse via SVD with a tolerance to zero-out small singular values.
# Input: X (m*n numeric matrix), tol (numeric threshold relative to max singular value).
# Output: n*m numeric matrix equal to X^+ (the pseudoinverse).
my_ginv <- function(X, tol = 1e-8) {
  s <- svd(X)
  d <- s$d
  d[ d < tol * max(d) ] <- 0  # zero out small singular values
  d_inv <- 1 / d
  d_inv[!is.finite(d_inv)] <- 0
  s$v %*% diag(d_inv, length(d_inv)) %*% t(s$u)
}


#################################################### This section simulates binary outcome data ###################################################

#Simulates grouped binary-response data under a logistic model with group-sparse coefficients.
#Input: n (samples), group.index (length p), s (group sparsity), beta.pre (optional k-vector).
#Output: list(X (n*p), Y (n), true.beta (p), A (p*k membership)).

simulate_data_binary=function(n,group.index,s=0,beta.pre=NULL){
  k=length(unique(group.index))
  p=length(group.index)
  A=matrix(0,nrow = p,k)
  for (i in 1:k) {
    A[,i]=(group.index==i)
  }
  if(is.null(beta.pre)){
    beta0=runif(k,1.5,2.5)
  }else{
    beta0=beta.pre
  }
  nonzero.ind=sample(1:k,ceiling(k*(1-s)))
  nonzero=rep(0,k)
  nonzero[nonzero.ind]=rep(c(-1,1),k%/%2)[1:length(nonzero.ind)]
  beta0=beta0*nonzero
  beta=A%*%beta0
  X=matrix(rpois(n*p,0.02),nrow = n,ncol = p)
  #X=matrix(rnorm(n*p,0.2,1),nrow = n,ncol = p)
  #X=t(apply(X, 1, function(x) x/sum(x)))
  means=crossprod(t(X),beta)
  ps=1/(1+exp(-means))
  Y=rbinom(n=n,size=1,prob = ps )
  #print(sigma)
  #sigma=1
  #X=X/sqrt(eigen(crossprod(X))$values[1])
  return(list(X=X,Y=Y,true.beta=beta,A=A))
}

############################# This section contains the simulation functions for 6.1.1 ####################################

# Finds all latent nodes whose subtree’s leaves EXACTLY match each user-specified leaf group.
# Input: tree_df (data.frame with columns node/parent/weight, etc.), leaf_group_vector (length = #leaves in root order from gather_leaves()).
# Output: named list mapping each group label to an integer vector of matching node ids (possibly empty).
find_latent_nodes_with_exact_leaf_groups <- function(tree_df, leaf_group_vector) {
  # ------------------------------------------------------------------
  # 1) Identify the root node and gather *all* leaves in the tree
  # ------------------------------------------------------------------
  root_node <- tree_df$node[is.na(tree_df$parent)]
  all_leaves <- gather_leaves(root_node, tree_df) 
  # all_leaves is sorted as per gather_leaves()
  
  # ------------------------------------------------------------------
  # 2) Match each leaf to the group index (leaf_group_vector)
  #    We assume leaf_group_vector is in the same order as 'all_leaves'.
  # ------------------------------------------------------------------
  leaf_groups_df <- data.frame(
    leaf  = all_leaves,
    group = leaf_group_vector
  )
  
  # Create a named list of leaf sets, one set per unique group index:
  # e.g., if group = c(1,1,2,2,3) you get:
  #   group_sets[["1"]] = c(leaf1, leaf2)
  #   group_sets[["2"]] = c(leaf3, leaf4)
  #   group_sets[["3"]] = c(leaf5)
  group_sets <- lapply(split(leaf_groups_df$leaf, leaf_groups_df$group), unique)
  
  # ------------------------------------------------------------------
  # 3) For each non-leaf node, retrieve the leaves in its subtree
  #    (the function gather_leaf_nodes_per_non_leaf() returns a data
  #     frame that has $node, $depth, and a list-column of $leaves, etc.)
  # ------------------------------------------------------------------
  gather_info     <- gather_leaf_nodes_per_non_leaf(tree_df)
  non_leaf_nodes  <- gather_info$df_result
  
  # non_leaf_nodes has columns:
  #   node, depth, leaves, latent_children, weight
  # where `leaves` is a list of leaf IDs for that node's subtree.
  
  # ------------------------------------------------------------------
  # 4) Check whether a node's subtree leaves exactly match a group set
  # ------------------------------------------------------------------
  same_set <- function(a, b) {
    if (length(a) != length(b)) return(FALSE)
    return(all(sort(a) == sort(b)))
  }
  
  # We'll create a named list: for each group, which nodes in the tree
  # have exactly that set of leaves?
  matching_latent_nodes <- vector("list", length(group_sets))
  names(matching_latent_nodes) <- names(group_sets)
  
  # For each group in group_sets:
  for (g in names(group_sets)) {
    group_leaves <- group_sets[[g]]
    
    # among all non-leaf nodes, find which have 'leaves' == group_leaves
    matched_nodes <- non_leaf_nodes$node[
      sapply(non_leaf_nodes$leaves, function(node_leaf_vec) {
        same_set(node_leaf_vec, group_leaves)
      })
    ]
    
    matching_latent_nodes[[g]] <- matched_nodes
  }
  
  # ------------------------------------------------------------------
  # 5) Return a list of latent nodes keyed by group index
  # ------------------------------------------------------------------
  return(matching_latent_nodes)
}

# Returns the full set of ancestor (parent, grandparent, …) node ids for a given set of nodes.
# Input: tree_df (data.frame with columns node and parent; rows indexed by node id), vec (integer vector of node ids).
# Output: integer vector of unique ancestor node ids (excluding NA/root), possibly empty.
find_all_parents=function(tree_df,vec){
  incre=T
  parent_vec=c()
  for(i in 1:length(vec)){
    parent_vec=c(parent_vec,tree_df$parent[vec[i]])
  }
  parent_vec=unique(parent_vec)
  ll=length(parent_vec)
  while(incre){
    for(i in 1:ll){
      parent_vec=c(parent_vec,tree_df$parent[parent_vec[i]])
    }
    parent_vec=unique(parent_vec)
    if(length(parent_vec)==ll) incre=F
    ll=length(parent_vec)
  }
  return(parent_vec[!is.na(parent_vec)])
}


# Runs a replicated simulation using a fixed tree and grouping: fits RARE, the tree-guided method, OLS, Ridge, and their oracle group versions; reports test MSEs and clustering agreement.
# Input: n (train size), n0 (test size), tree_df (node/parent/weight data.frame for the tree), tree_group (length-p group labels), n1 (validation size), s (sparsity), reps (#replicates), ridge.param, thresh, INratio (signal-to-noise control).
# Output: list of length-'reps' vectors: our.minloss, rare.minloss, our.rand, rare.rand, our.minloss.ideal, rare.minloss.ideal, our.rand.ideal, rare.rand.ideal, oracle.ls.loss, ls.loss, ridge.loss, oracle.ridge.loss.
simulate_error_given_tree=function(n,n0, trees,tree_df,tree_group,n1=n,s=0,reps=50,ridge.param=0,thresh=1e-5,INratio=5){
  rare.minloss=numeric(reps)
  our.minloss=numeric(reps)
  rare.rand=numeric(reps)
  our.rand=numeric(reps)
  rare.minloss.ideal=numeric(reps)
  our.minloss.ideal=numeric(reps)
  rare.rand.ideal=numeric(reps)
  our.rand.ideal=numeric(reps)
  oracle.ls.loss=numeric(reps)
  ls.loss=numeric(reps)
  ridge.loss         = numeric(reps)
  oracle.ridge.loss  = numeric(reps)
  for (i in 1:reps) {
    cat(i,"/",reps,'\n')
    data=simulate_data(n+n0+n1,tree_group,ratio = INratio,scale.X=TRUE)
    train_index=sample(1:(n+n0+n1),n)
    valid_index=sample((1:(n+n0+n1))[- train_index],n1)
    test_index=(1:(n+n0+n1))[- c(train_index,valid_index)]
    
    rare.result=rarefit(data$Y[train_index],data$X[train_index,],hc = trees$tree,alpha = 1,intercept = F)
    loss1=apply(data$X[valid_index,]%*%rare.result$beta[[1]],2,function(x) sum((x-data$Y[valid_index])^2))
    loss3=apply(data$X[valid_index,]%*%rare.result$beta[[1]],2,function(x) sum((x-data$X[valid_index,]%*%data$true.beta)^2))
    
    rare.beta=rare.result$beta[[1]][,which.min(loss1)]
    rare.beta.ideal=rare.result$beta[[1]][,which.min(loss3)]
    rare.minloss[i]=sum((data$Y[test_index]-data$X[test_index,]%*%rare.beta)^2)/n0
    rare.minloss.ideal[i]=sum((data$Y[test_index]-data$X[test_index,]%*%rare.beta.ideal)^2)/n0
    our.result=grid.simple_linear(Y = data$Y[train_index],X=data$X[train_index,],tree_df,true_beta= data$true.beta,ridge.param,thresh = thresh)
    loss2=apply(data$X[valid_index,]%*%our.result$beta,2,function(x) sum((x-data$Y[valid_index])^2))
    loss4=apply(data$X[valid_index,]%*%our.result$beta,2,function(x) sum((x-data$X[valid_index,]%*%data$true.beta)^2))
    our.beta=our.result$beta[,which.min(loss2)]
    our.beta.ideal=our.result$beta[,which.min(loss4)]
    our.minloss[i]=sum((data$Y[test_index]-data$X[test_index,]%*%our.beta)^2)/n0
    our.minloss.ideal[i]=sum((data$Y[test_index]-data$X[test_index,]%*%our.beta.ideal)^2)/n0
    rare.rand[i]=adjustedRandIndex(as.numeric(as.factor(data$true.beta)),as.numeric(as.factor(rare.beta)))
    rare.rand.ideal[i]=adjustedRandIndex(as.numeric(as.factor(data$true.beta)),as.numeric(as.factor(rare.beta.ideal)))
    our.rand[i]=adjustedRandIndex(as.numeric(as.factor(data$true.beta)),as.numeric(as.factor(our.beta)))
    our.rand.ideal[i]=adjustedRandIndex(as.numeric(as.factor(data$true.beta)),as.numeric(as.factor(our.beta.ideal)))
    
    
    ls.X=data$X[train_index,]
    ls.Y=data$Y[train_index]
    ls.ginv=my_ginv(crossprod(ls.X))
    ls.coef=crossprod(ls.ginv,crossprod(ls.X,ls.Y))
    ls.loss[i]=sum((data$Y[test_index]-data$X[test_index,]%*%ls.coef)^2)/n0
    
    Q=make_new_X(data$X,trees$group)
    ols.X=Q$X1[train_index,]
    ols.ginv=t(my_ginv(crossprod(ols.X)))
    ols.coef=crossprod(ols.ginv,crossprod(ols.X,ls.Y))
    oracle.ls.loss[i]=sum((data$Y[test_index]-Q$X1[test_index,]%*%ols.coef)^2)/n0
    
    X_train = data$X[train_index,]
    Y_train = data$Y[train_index]
    X_valid = data$X[valid_index,]
    Y_valid = data$Y[valid_index]
    X_test  = data$X[test_index,]
    Y_test  = data$Y[test_index]
    
    # alpha=0 => Ridge, intercept=FALSE => no intercept in the model
    fit_ridge = glmnet(X_train, Y_train, alpha = 0, intercept = FALSE)
    
    # The sequence of lambdas in fit_ridge
    lam_seq   = fit_ridge$lambda
    
    # Predictions on validation set for each lambda
    preds_valid_ridge = predict(fit_ridge, newx = X_valid, s = lam_seq)
    # Compute validation MSE for each column (each lambda)
    val_mse_ridge = colMeans((Y_valid - preds_valid_ridge)^2)
    
    best_lam_std  = lam_seq[which.min(val_mse_ridge)]
    
    # Extract coefficients at best lambda (no intercept in the model)
    # 'coef(...)' returns a vector of length p+1 (the first is intercept).
    # But intercept=FALSE means that will be zero, so we can exclude it.
    beta_ridge_best_full = as.numeric(coef(fit_ridge, s = best_lam_std))
    beta_ridge_best      = beta_ridge_best_full[-1]  # omit intercept slot
    
    # Evaluate on test
    preds_ridge_test = X_test %*% beta_ridge_best  # no intercept
    ridge.loss[i] = mean((Y_test - preds_ridge_test)^2)
    
    ##-----------------------------------------------------------------
    ## 6. NEW: Oracle Ridge using glmnet on grouped design Q$X1
    ##         with no intercept
    ##-----------------------------------------------------------------
    or.X_train = Q$X1[train_index, ]
    or.X_valid = Q$X1[valid_index, ]
    or.X_test  = Q$X1[test_index, ]
    
    fit_oracle_ridge = glmnet(or.X_train, Y_train, alpha = 0, intercept = FALSE)
    lam_seq_or       = fit_oracle_ridge$lambda
    
    preds_valid_or = predict(fit_oracle_ridge, newx = or.X_valid, s = lam_seq_or)
    val_mse_oracle_ridge = colMeans((Y_valid - preds_valid_or)^2)
    
    best_lam_oracle = lam_seq_or[which.min(val_mse_oracle_ridge)]
    
    # Extract group-level coefficients at best lambda
    beta_group_best_full = as.numeric(coef(fit_oracle_ridge, s = best_lam_oracle))
    beta_group_best      = beta_group_best_full[-1]  # omit the intercept slot
    
    # Map from group-level to full p-dim coefficients in feature space
    beta_or_feature = Q$Q %*% beta_group_best
    
    # Evaluate on test set (no intercept)
    preds_test_or = X_test %*% beta_or_feature
    oracle.ridge.loss[i] = mean((Y_test - preds_test_or)^2)
    
    
  }
  return(list(our.minloss=our.minloss,rare.minloss=rare.minloss,our.rand=our.rand,rare.rand=rare.rand,our.minloss.ideal=our.minloss.ideal,rare.minloss.ideal=rare.minloss.ideal,our.rand.ideal=our.rand.ideal,rare.rand.ideal=rare.rand.ideal,oracle.ls.loss=oracle.ls.loss,ls.loss=ls.loss,ridge.loss = ridge.loss,oracle.ridge.loss = oracle.ridge.loss))#,oracle.ridge.loss=oracle.ridge.loss))
}             


############################# This section contains the simulation functions for 6.1.2 ####################################


# Generates a block-structured distance matrix for k equal-sized subtrees (p0 leaves each) and builds an hclust tree.
# Input: k (number of subtrees/groups), p0 (leaves per subtree).
# Output: list(tree = hclust object from the combined distances, group = length-(k*p0) integer vector of group labels 1...k).
simulate_combined_tree <- function(k, p0) {
  # Step 1: Initialize variables
  total_leaves <- k * p0
  group_labels <- c()
  combined_dist <- matrix(0, nrow = total_leaves, ncol = total_leaves)
  
  # Step 2: Generate k subtrees
  start_index <- 1
  for (i in 1:k) {
    # Generate latent values for the current subtree
    latent_values <- rnorm(p0, mean = i, sd = 0.05)
    
    # Compute distances for the current subtree
    subtree_dist <- as.matrix(dist(latent_values))
    
    # Place subtree distances into the combined matrix
    end_index <- start_index + p0 - 1
    combined_dist[start_index:end_index, start_index:end_index] <- subtree_dist
    
    # Assign group labels
    group_labels <- c(group_labels, rep(i, p0))
    
    start_index <- end_index + 1
  }
  
  # Step 3: Add fixed distances between subtrees
  combined_dist[combined_dist == 0] <- Inf  # Temporarily mark empty cells
  block_indices <- seq(1, total_leaves, by = p0)
  
  for (i in 1:(k - 1)) {
    start1 <- block_indices[i]
    end1 <- start1 + p0 - 1
    start2 <- block_indices[i + 1]
    end2 <- start2 + p0 - 1
    
    combined_dist[start1:end1, start2:end2] <- 10  # Fixed distance between subtrees
    combined_dist[start2:end2, start1:end1] <- 10  # Symmetric assignment
  }
  
  # Step 4: Replace remaining Inf values with large finite distances
  combined_dist[is.infinite(combined_dist)] <- max(combined_dist[is.finite(combined_dist)]) * 2
  
  # Step 5: Convert to a valid distance object
  combined_dist <- as.dist(combined_dist)
  
  # Step 6: Perform hierarchical clustering
  combined_tree <- hclust(combined_dist, method = "average")
  
  # Return the hierarchical tree and group labels
  return(list(tree = combined_tree, group = group_labels))
}

# Repeatedly simulates data under a fixed/provided tree, fits multiple methods (RARE, ours, OLS/Ridge + oracle variants), and evaluates test MSE and ARI.
# Input: n (train size), n0 (test size), p0 (leaves per group), k (number of groups), n1 (validation size, default n), s (sparsity), reps, ridge.param, thresh, beta.specify (optional group-level β), tree.specify (optional prebuilt list(tree, group)), INratio (SNR), weight.order.
# Output: list of length-reps vectors: our.minloss, rare.minloss, our.rand, rare.rand, our.minloss.ideal, rare.minloss.ideal, our.rand.ideal, rare.rand.ideal, oracle.ls.loss, ls.loss, ridge.loss, oracle.ridge.loss.
simulate_error_tree1_valid=function(n,n0,p0,k,n1=n,s=0,reps=50,ridge.param=0,thresh=1e-5,beta.specify=NULL,tree.specify=NULL,INratio=5,weight.order=-1,fix_sigma=NULL){
  rare.minloss=numeric(reps)
  our.minloss=numeric(reps)
  rare.rand=numeric(reps)
  our.rand=numeric(reps)
  rare.minloss.ideal=numeric(reps)
  our.minloss.ideal=numeric(reps)
  rare.rand.ideal=numeric(reps)
  our.rand.ideal=numeric(reps)
  oracle.ls.loss=numeric(reps)
  ls.loss=numeric(reps)
  ridge.loss         = numeric(reps)
  oracle.ridge.loss  = numeric(reps)
  if(is.null(tree.specify)) trees=simulate_combined_tree(k,p0)
  else trees=tree.specify
  for (i in 1:reps) {
    cat(i,"/",reps,'\n')
    tree_df=hclust_to_df(trees$tree)
    data=simulate_data(n+n0+n1,trees$group,beta.pre = beta.specify,ratio=INratio,scale.X=TRUE,fix_sigma=fix_sigma)
    train_index=sample(1:(n+n0+n1),n)
    valid_index=sample((1:(n+n0+n1))[- train_index],n1)
    test_index=(1:(n+n0+n1))[- c(train_index,valid_index)]
    
    rare.result=rarefit(data$Y[train_index],data$X[train_index,],hc = trees$tree,alpha = 1,intercept = F)
    loss1=apply(data$X[valid_index,]%*%rare.result$beta[[1]],2,function(x) sum((x-data$Y[valid_index])^2))
    loss3=apply(data$X[valid_index,]%*%rare.result$beta[[1]],2,function(x) sum((x-data$X[valid_index,]%*%data$true.beta)^2))
    
    rare.beta=rare.result$beta[[1]][,which.min(loss1)]
    rare.beta.ideal=rare.result$beta[[1]][,which.min(loss3)]
    rare.minloss[i]=sum((data$Y[test_index]-data$X[test_index,]%*%rare.beta)^2)/n0
    rare.minloss.ideal[i]=sum((data$Y[test_index]-data$X[test_index,]%*%rare.beta.ideal)^2)/n0
    our.result=grid.simple_linear(Y = data$Y[train_index],X=data$X[train_index,],tree_df,true_beta= data$true.beta,ridge.param,thresh = thresh)
    loss2=apply(data$X[valid_index,]%*%our.result$beta,2,function(x) sum((x-data$Y[valid_index])^2))
    loss4=apply(data$X[valid_index,]%*%our.result$beta,2,function(x) sum((x-data$X[valid_index,]%*%data$true.beta)^2))
    our.beta=our.result$beta[,which.min(loss2)]
    our.beta.ideal=our.result$beta[,which.min(loss4)]
    our.minloss[i]=sum((data$Y[test_index]-data$X[test_index,]%*%our.beta)^2)/n0
    our.minloss.ideal[i]=sum((data$Y[test_index]-data$X[test_index,]%*%our.beta.ideal)^2)/n0
    rare.rand[i]=adjustedRandIndex(as.numeric(as.factor(data$true.beta)),as.numeric(as.factor(rare.beta)))
    rare.rand.ideal[i]=adjustedRandIndex(as.numeric(as.factor(data$true.beta)),as.numeric(as.factor(rare.beta.ideal)))
    our.rand[i]=adjustedRandIndex(as.numeric(as.factor(data$true.beta)),as.numeric(as.factor(our.beta)))
    our.rand.ideal[i]=adjustedRandIndex(as.numeric(as.factor(data$true.beta)),as.numeric(as.factor(our.beta.ideal)))
    
    
    ls.X=data$X[train_index,]
    ls.Y=data$Y[train_index]
    ls.ginv=my_ginv(crossprod(ls.X))
    ls.coef=crossprod(ls.ginv,crossprod(ls.X,ls.Y))
    ls.loss[i]=sum((data$Y[test_index]-data$X[test_index,]%*%ls.coef)^2)/n0
    
    Q=make_new_X(data$X,trees$group)
    ols.X=Q$X1[train_index,]
    ols.ginv=t(my_ginv(crossprod(ols.X)))
    ols.coef=crossprod(ols.ginv,crossprod(ols.X,ls.Y))
    oracle.ls.loss[i]=sum((data$Y[test_index]-Q$X1[test_index,]%*%ols.coef)^2)/n0
    
    X_train = data$X[train_index,]
    Y_train = data$Y[train_index]
    X_valid = data$X[valid_index,]
    Y_valid = data$Y[valid_index]
    X_test  = data$X[test_index,]
    Y_test  = data$Y[test_index]
    
    # alpha=0 => Ridge, intercept=FALSE => no intercept in the model
    fit_ridge = glmnet(X_train, Y_train, alpha = 0, intercept = FALSE)
    
    # The sequence of lambdas in fit_ridge
    lam_seq   = fit_ridge$lambda
    
    # Predictions on validation set for each lambda
    preds_valid_ridge = predict(fit_ridge, newx = X_valid, s = lam_seq)
    # Compute validation MSE for each column (each lambda)
    val_mse_ridge = colMeans((Y_valid - preds_valid_ridge)^2)
    
    best_lam_std  = lam_seq[which.min(val_mse_ridge)]
    
    # Extract coefficients at best lambda (no intercept in the model)
    # 'coef(...)' returns a vector of length p+1 (the first is intercept).
    # But intercept=FALSE means that will be zero, so we can exclude it.
    beta_ridge_best_full = as.numeric(coef(fit_ridge, s = best_lam_std))
    beta_ridge_best      = beta_ridge_best_full[-1]  # omit intercept slot
    
    # Evaluate on test
    preds_ridge_test = X_test %*% beta_ridge_best  # no intercept
    ridge.loss[i] = mean((Y_test - preds_ridge_test)^2)
    
    ##-----------------------------------------------------------------
    ## 6. NEW: Oracle Ridge using glmnet on grouped design Q$X1
    ##         with no intercept
    ##-----------------------------------------------------------------
    or.X_train = Q$X1[train_index, ]
    or.X_valid = Q$X1[valid_index, ]
    or.X_test  = Q$X1[test_index, ]
    
    fit_oracle_ridge = glmnet(or.X_train, Y_train, alpha = 0, intercept = FALSE)
    lam_seq_or       = fit_oracle_ridge$lambda
    
    preds_valid_or = predict(fit_oracle_ridge, newx = or.X_valid, s = lam_seq_or)
    val_mse_oracle_ridge = colMeans((Y_valid - preds_valid_or)^2)
    
    best_lam_oracle = lam_seq_or[which.min(val_mse_oracle_ridge)]
    
    # Extract group-level coefficients at best lambda
    beta_group_best_full = as.numeric(coef(fit_oracle_ridge, s = best_lam_oracle))
    beta_group_best      = beta_group_best_full[-1]  # omit the intercept slot
    
    # Map from group-level to full p-dim coefficients in feature space
    beta_or_feature = Q$Q %*% beta_group_best
    
    # Evaluate on test set (no intercept)
    preds_test_or = X_test %*% beta_or_feature
    oracle.ridge.loss[i] = mean((Y_test - preds_test_or)^2)
    
    
  }
  return(list(our.minloss=our.minloss,rare.minloss=rare.minloss,our.rand=our.rand,rare.rand=rare.rand,our.minloss.ideal=our.minloss.ideal,rare.minloss.ideal=rare.minloss.ideal,our.rand.ideal=our.rand.ideal,rare.rand.ideal=rare.rand.ideal,oracle.ls.loss=oracle.ls.loss,ls.loss=ls.loss,ridge.loss = ridge.loss,oracle.ridge.loss = oracle.ridge.loss))#,oracle.ridge.loss=oracle.ridge.loss))
}      


##################### This section contains simulation functions for 6.2_1 #############################

# Runs repeated simulations comparing RARE, our tree-guided method, OLS/Ridge, and oracle variants; evaluates test MSE and clustering agreement.
# Input: n (train size), n0 (test size), p (features), k (#groups), n1 (validation size), s (sparsity), reps, ridge.param, thresh, INratio (SNR), weight.order.
# Output: list of length-reps vectors: our.minloss, rare.minloss, our.rand, rare.rand, our.minloss.ideal, rare.minloss.ideal, our.rand.ideal, rare.rand.ideal, oracle.ls.loss, ls.loss, ridge.loss, oracle.ridge.loss.
simulate_error_valid=function(n,n0,p,k,n1=n,s=0,reps=50,ridge.param=0,thresh=1e-5,INratio=5,weight.order=-1){
  rare.minloss=numeric(reps)
  our.minloss=numeric(reps)
  rare.rand=numeric(reps)
  our.rand=numeric(reps)
  rare.minloss.ideal=numeric(reps)
  our.minloss.ideal=numeric(reps)
  rare.rand.ideal=numeric(reps)
  our.rand.ideal=numeric(reps)
  oracle.ls.loss=numeric(reps)
  ls.loss=numeric(reps)
  ridge.loss         = numeric(reps)
  oracle.ridge.loss  = numeric(reps)
  for (i in 1:reps) {
    cat(i,"/",reps,'\n')
    trees=simulate_tree(k,p)
    tree_df=hclust_to_df(trees$tree,weight.order=weight.order)
    data=simulate_data(n+n0+n1,trees$group,ratio = INratio)
    train_index=sample(1:(n+n0+n1),n)
    valid_index=sample((1:(n+n0+n1))[- train_index],n1)
    test_index=(1:(n+n0+n1))[- c(train_index,valid_index)]
    
    rare.result=rarefit(data$Y[train_index],data$X[train_index,],hc = trees$tree,alpha = 1,intercept = F)
    loss1=apply(data$X[valid_index,]%*%rare.result$beta[[1]],2,function(x) sum((x-data$Y[valid_index])^2))
    loss3=apply(data$X[valid_index,]%*%rare.result$beta[[1]],2,function(x) sum((x-data$X[valid_index,]%*%data$true.beta)^2))
    
    rare.beta=rare.result$beta[[1]][,which.min(loss1)]
    rare.beta.ideal=rare.result$beta[[1]][,which.min(loss3)]
    rare.minloss[i]=sum((data$Y[test_index]-data$X[test_index,]%*%rare.beta)^2)/n0
    rare.minloss.ideal[i]=sum((data$Y[test_index]-data$X[test_index,]%*%rare.beta.ideal)^2)/n0
    our.result=grid.simple_linear(Y = data$Y[train_index],X=data$X[train_index,],tree_df,true_beta= data$true.beta,ridge.param,thresh = thresh)
    loss2=apply(data$X[valid_index,]%*%our.result$beta,2,function(x) sum((x-data$Y[valid_index])^2))
    loss4=apply(data$X[valid_index,]%*%our.result$beta,2,function(x) sum((x-data$X[valid_index,]%*%data$true.beta)^2))
    our.beta=our.result$beta[,which.min(loss2)]
    our.beta.ideal=our.result$beta[,which.min(loss4)]
    our.minloss[i]=sum((data$Y[test_index]-data$X[test_index,]%*%our.beta)^2)/n0
    our.minloss.ideal[i]=sum((data$Y[test_index]-data$X[test_index,]%*%our.beta.ideal)^2)/n0
    rare.rand[i]=adjustedRandIndex(as.numeric(as.factor(data$true.beta)),as.numeric(as.factor(rare.beta)))
    rare.rand.ideal[i]=adjustedRandIndex(as.numeric(as.factor(data$true.beta)),as.numeric(as.factor(rare.beta.ideal)))
    our.rand[i]=adjustedRandIndex(as.numeric(as.factor(data$true.beta)),as.numeric(as.factor(our.beta)))
    our.rand.ideal[i]=adjustedRandIndex(as.numeric(as.factor(data$true.beta)),as.numeric(as.factor(our.beta.ideal)))
    
    
    ls.X=data$X[train_index,]
    ls.Y=data$Y[train_index]
    ls.ginv=my_ginv(crossprod(ls.X))
    ls.coef=crossprod(ls.ginv,crossprod(ls.X,ls.Y))
    ls.loss[i]=sum((data$Y[test_index]-data$X[test_index,]%*%ls.coef)^2)/n0
    
    Q=make_new_X(data$X,trees$group)
    ols.X=Q$X1[train_index,]
    ols.ginv=t(my_ginv(crossprod(ols.X)))
    ols.coef=crossprod(ols.ginv,crossprod(ols.X,ls.Y))
    oracle.ls.loss[i]=sum((data$Y[test_index]-Q$X1[test_index,]%*%ols.coef)^2)/n0
    
    X_train = data$X[train_index,]
    Y_train = data$Y[train_index]
    X_valid = data$X[valid_index,]
    Y_valid = data$Y[valid_index]
    X_test  = data$X[test_index,]
    Y_test  = data$Y[test_index]
    
    # alpha=0 => Ridge, intercept=FALSE => no intercept in the model
    fit_ridge = glmnet(X_train, Y_train, alpha = 0, intercept = FALSE)
    
    # The sequence of lambdas in fit_ridge
    lam_seq   = fit_ridge$lambda
    
    # Predictions on validation set for each lambda
    preds_valid_ridge = predict(fit_ridge, newx = X_valid, s = lam_seq)
    # Compute validation MSE for each column (each lambda)
    val_mse_ridge = colMeans((Y_valid - preds_valid_ridge)^2)
    
    best_lam_std  = lam_seq[which.min(val_mse_ridge)]
    
    # Extract coefficients at best lambda (no intercept in the model)
    # 'coef(...)' returns a vector of length p+1 (the first is intercept).
    # But intercept=FALSE means that will be zero, so we can exclude it.
    beta_ridge_best_full = as.numeric(coef(fit_ridge, s = best_lam_std))
    beta_ridge_best      = beta_ridge_best_full[-1]  # omit intercept slot
    
    # Evaluate on test
    preds_ridge_test = X_test %*% beta_ridge_best  # no intercept
    ridge.loss[i] = mean((Y_test - preds_ridge_test)^2)
    
    ##-----------------------------------------------------------------
    ## 6. NEW: Oracle Ridge using glmnet on grouped design Q$X1
    ##         with no intercept
    ##-----------------------------------------------------------------
    or.X_train = Q$X1[train_index, ]
    or.X_valid = Q$X1[valid_index, ]
    or.X_test  = Q$X1[test_index, ]
    
    fit_oracle_ridge = glmnet(or.X_train, Y_train, alpha = 0, intercept = FALSE)
    lam_seq_or       = fit_oracle_ridge$lambda
    
    preds_valid_or = predict(fit_oracle_ridge, newx = or.X_valid, s = lam_seq_or)
    val_mse_oracle_ridge = colMeans((Y_valid - preds_valid_or)^2)
    
    best_lam_oracle = lam_seq_or[which.min(val_mse_oracle_ridge)]
    
    # Extract group-level coefficients at best lambda
    beta_group_best_full = as.numeric(coef(fit_oracle_ridge, s = best_lam_oracle))
    beta_group_best      = beta_group_best_full[-1]  # omit the intercept slot
    
    # Map from group-level to full p-dim coefficients in feature space
    beta_or_feature = Q$Q %*% beta_group_best
    
    # Evaluate on test set (no intercept)
    preds_test_or = X_test %*% beta_or_feature
    oracle.ridge.loss[i] = mean((Y_test - preds_test_or)^2)
    
    
  }
  return(list(our.minloss=our.minloss,rare.minloss=rare.minloss,our.rand=our.rand,rare.rand=rare.rand,our.minloss.ideal=our.minloss.ideal,rare.minloss.ideal=rare.minloss.ideal,our.rand.ideal=our.rand.ideal,rare.rand.ideal=rare.rand.ideal,oracle.ls.loss=oracle.ls.loss,ls.loss=ls.loss,ridge.loss = ridge.loss,oracle.ridge.loss = oracle.ridge.loss))#,oracle.ridge.loss=oracle.ridge.loss))
}      


########################################################## This section simulates Section 6.2_3 #############################################################

#Repeatedly simulates data, fits RARE and the tree-guided model (plus OLS/Ridge and oracle variants), and evaluates test losses and ARI.
#Input: n (train), n0 (test), p, k (#groups), n1 (valid), s, reps, ridge.param, thresh, weight.order.
#Output: list of length-reps vectors: our.minloss, rare.minloss, our.rand, rare.rand, our.minloss.ideal, rare.minloss.ideal, our.rand.ideal, rare.rand.ideal, oracle.ls.loss, ls.loss, ridge.loss, oracle.ridge.loss.
simulate_error_valid.logistic=function(n,n0,p,k,n1=n,s=0,reps=50,ridge.param=0,thresh=1e-5,weight.order=-1){
  n1=n
  rare.minloss=numeric(reps)
  our.minloss=numeric(reps)
  rare.rand=numeric(reps)
  our.rand=numeric(reps)
  rare.minloss.ideal=numeric(reps)
  our.minloss.ideal=numeric(reps)
  rare.rand.ideal=numeric(reps)
  our.rand.ideal=numeric(reps)
  oracle.ls.loss=numeric(reps)
  ls.loss=numeric(reps)
  ridge.loss         = numeric(reps)
  oracle.ridge.loss  = numeric(reps)
  for (i in 1:reps) {
    cat(i,"/",reps,'\n')
    trees=simulate_tree(k,p)
    tree_df=hclust_to_df(trees$tree,weight.order=weight.order)
    data=simulate_data_binary(n = n+n0+n1,group.index = trees$group)
    train_index=sample(1:(n+n0+n1),n)
    valid_index=sample((1:(n+n0+n1))[- train_index],n1)
    test_index=(1:(n+n0+n1))[- c(train_index,valid_index)]
    
    tree_df = hclust_to_df(trees$tree, weight.order=weight.order)
    rare.result=rarefit.logistic(data$Y[train_index],data$X[train_index,],tree_df = tree_df,intercept = F)
    loss1=negtv_lglkh(data$Y[valid_index],data$X[valid_index,],rare.result$beta)
    loss3=negtv_lglkh(1/(1+exp(-crossprod(t(data$X[valid_index,]),data$true.beta))),data$X[valid_index,],rare.result$beta)
    
    rare.beta=rare.result$beta[,which.min(loss1)]
    rare.beta.ideal=rare.result$beta[,which.min(loss3)]
    rare.minloss[i]=negtv_lglkh(data$Y[test_index],data$X[test_index,],rare.beta)
    rare.minloss.ideal[i]=negtv_lglkh(data$Y[test_index],data$X[test_index,],rare.beta.ideal)
    
    our.result=grid.logistic(Y = data$Y[train_index],X=data$X[train_index,],tree_df,true_beta= data$true.beta,ridge.param,thresh = thresh)
    loss2=negtv_lglkh(data$Y[valid_index],data$X[valid_index,],our.result$beta)
    loss4=negtv_lglkh(1/(1+exp(-crossprod(t(data$X[valid_index,]),data$true.beta))),data$X[valid_index,],our.result$beta)
    our.beta=our.result$beta[,which.min(loss2)]
    our.beta.ideal=our.result$beta[,which.min(loss4)]
    our.minloss[i]=negtv_lglkh(data$Y[test_index],data$X[test_index,],our.beta)
    our.minloss.ideal[i]=negtv_lglkh(data$Y[test_index],data$X[test_index,],our.beta.ideal)
    
    rare.rand[i]=adjustedRandIndex(as.numeric(as.factor(data$true.beta)),as.numeric(as.factor(rare.beta)))
    rare.rand.ideal[i]=adjustedRandIndex(as.numeric(as.factor(data$true.beta)),as.numeric(as.factor(rare.beta.ideal)))
    our.rand[i]=adjustedRandIndex(as.numeric(as.factor(data$true.beta)),as.numeric(as.factor(our.beta)))
    our.rand.ideal[i]=adjustedRandIndex(as.numeric(as.factor(data$true.beta)),as.numeric(as.factor(our.beta.ideal)))
    
    
    ls.X=data$X[train_index,]
    ls.Y=data$Y[train_index]
    ls.res=glmnet(y=ls.Y,x=ls.X,lambda = 0,intercept = FALSE,nlambda = 50)
    ls.coef=as.numeric(ls.res$beta)
    
    ls.loss[i]=negtv_lglkh(data$Y[test_index],data$X[test_index,],ls.coef)
    
    Q=make_new_X(data$X,trees$group)
    ols.X=Q$X1[train_index,]
    ols.res=glmnet(y=ls.Y,x=ols.X,lambda = 0,intercept = FALSE,nlambda = 50)
    ols.coef=as.numeric(ols.res$beta)
    
    oracle.ls.loss[i]=negtv_lglkh(data$Y[test_index],Q$X1[test_index,],ols.coef)
    
    
    X_train = data$X[train_index,]
    Y_train = data$Y[train_index]
    X_valid = data$X[valid_index,]
    Y_valid = data$Y[valid_index]
    X_test  = data$X[test_index,]
    Y_test  = data$Y[test_index]
    
    # alpha=0 => Ridge, intercept=FALSE => no intercept in the model
    fit_ridge = glmnet(X_train, Y_train, alpha = 0, intercept = FALSE,family = "binomial",nlambda = 50)
    
    
    val_mse_ridge = negtv_lglkh(Y_valid,X_valid,as.matrix(fit_ridge$beta))
    
    ridge.beta=fit_ridge$beta[,which.min(val_mse_ridge)]
    
    ridge.loss[i] = negtv_lglkh(Y_test,X_test,ridge.beta)
    
    ##-----------------------------------------------------------------
    ## 6. NEW: Oracle Ridge using glmnet on grouped design Q$X1
    ##         with no intercept
    ##-----------------------------------------------------------------
    or.X_train = Q$X1[train_index, ]
    or.X_valid = Q$X1[valid_index, ]
    or.X_test  = Q$X1[test_index, ]
    
    fit_oracle_ridge = glmnet(or.X_train, Y_train, alpha = 0, intercept = FALSE,family = "binomial",nlambda = 50)
    val_mse_oracle_ridge = negtv_lglkh(Y_valid,or.X_valid,as.matrix(fit_oracle_ridge$beta))
    or.ridge.beta=fit_oracle_ridge$beta[,which.min(val_mse_oracle_ridge)]
    
    oracle.ridge.loss[i] = negtv_lglkh(Y_test,or.X_test,or.ridge.beta)
    
  
    #bet=Q$Q%*%(solve(t(Q$X1)%*%Q$X1+ridge.param*diag(ncol(Q$Q)))%*%t(Q$X1)%*%data$Y)
    #oracle.ridge.loss[i]=sum((bet-data$true.beta)^2)/p
  }
  return(list(our.minloss=our.minloss,rare.minloss=rare.minloss,our.rand=our.rand,rare.rand=rare.rand,our.minloss.ideal=our.minloss.ideal,rare.minloss.ideal=rare.minloss.ideal,our.rand.ideal=our.rand.ideal,rare.rand.ideal=rare.rand.ideal,oracle.ls.loss=oracle.ls.loss,ls.loss=ls.loss,ridge.loss = ridge.loss,oracle.ridge.loss = oracle.ridge.loss))#,oracle.ridge.loss=oracle.ridge.loss))
}
