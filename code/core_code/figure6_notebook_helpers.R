# Original notebook helpers used to reproduce Figure 6.
# Function definitions are copied verbatim, retaining the last definition
# encountered when the two original R files are sourced in notebook order.
# Original source files and SHA-256 hashes:
#   functions_v2.R: a655e53b48de51ce255fff6bd80fb4789ea37398dd83bb90227572ba0cfdd9dd
#   logistic_functions_v2.R: 6c484f9b5f20d7116edf18b178caafa57f6410b4a566f2aecf162234bfa81cdf
#
# The entry script loads glmnet and compiles figure6_notebook_algos.cpp
# to provide prox_tree before calling these functions. The C++ file is a
# byte-for-byte copy of direct_feature_aggregation/code/algos.cpp, the
# original proximal implementation present in the repository:
#   algos.cpp SHA-256: 0ce481edbd3ddcb775a9a251758b06fd66de5ce2822aa836fd24d000b865c581
# No package loading, compilation, fitting, or path changes occur here.

# Source: functions_v2.R:23
calculate_depth <- function(node, df) {
  depth <- 0
  while (!is.na(node)) {
    node <- df$parent[df$node == node]
    depth <- depth + 1
  }
  return(depth - 1)
}

# Source: functions_v2.R:32
detect_non_leaf_nodes <- function(df) {
  # Identify all nodes that are not leaf nodes
  non_leaf_nodes <- unique(df$parent[!is.na(df$parent)])
  
  # Calculate the depth of each non-leaf node
  non_leaf_depths <- sapply(non_leaf_nodes, calculate_depth, df = df)
  
  # Combine non-leaf nodes and their depths into a data frame
  non_leaf_df <- data.frame(node = non_leaf_nodes, depth = non_leaf_depths)
  
  return(non_leaf_df)
}

# Source: functions_v2.R:45
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

# Source: functions_v2.R:68
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

# Source: functions_v2.R:80
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

# Source: functions_v2.R:107
gather_leaf_nodes_per_non_leaf <- function(df) {
  # Detect non-leaf nodes and their depths
  non_leaf_nodes_df <- detect_non_leaf_nodes(df)
  
  # Gather leaf nodes, latent children, and weights for each non-leaf node
  non_leaf_nodes_df$leaves <- lapply(non_leaf_nodes_df$node, gather_leaves, df = df)
  non_leaf_nodes_df$latent_children <- lapply(non_leaf_nodes_df$node, gather_direct_latent_nodes, df = df)
  non_leaf_nodes_df$weight <- sapply(non_leaf_nodes_df$node, function(node) df$weight[df$node == node])
  
  result <- make_list(non_leaf_nodes_df)
  return(list(list_result = result, df_result = non_leaf_nodes_df))
}

# Source: functions_v2.R:200
find_p=function(tree_df){
    roots=tree_df$node[is.na(tree_df$parent)]
    p=0
    for(i in 1:length(roots)){
        p=p+length(gather_leaves(roots[i],tree_df))
    }
    return(p)
}

# Source: functions_v2.R:208
find_coarest=function (df, result) 
{
    p = find_p(df)
    all_cover = F
    max_depth = max(result$depth)
    coarest_set = result[result$depth == 0, ]
    if (max(result[result$depth == 0, ]$weight) != 0) {
        indiv.nodes = (1:p)[-sort(unlist((coarest_set$leaves)))]
        if (length(indiv.nodes) != 0) {
            for (i in 1:length(indiv.nodes)) {
                new_df = data.frame(node = indiv.nodes[i], depth = 1, 
                  leaves = list(indiv.nodes[i]), latent_children = NA, 
                  weight = 0)
                names(new_df) = names(coarest_set)
                coarest_set = rbind(coarest_set, new_df)
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
            indiv.nodes = (1:p)[-sort(unlist((coarest_set$leaves)))]
            if (length(indiv.nodes) != 0) {
                for (i in 1:length(indiv.nodes)) {
                  new_df = data.frame(node = indiv.nodes[i], 
                    depth = 1, leaves = list(indiv.nodes[i]), 
                    latent_children = NA, weight = 0)
                  names(new_df) = names(coarest_set)
                  coarest_set = rbind(coarest_set, new_df)
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
    indiv.nodes = (1:p)[-sort(unlist((coarest_set$leaves)))]
    if (length(indiv.nodes) != 0) {
        for (i in 1:length(indiv.nodes)) {
            new_df = data.frame(node = indiv.nodes[i], depth = 1, 
                leaves = list(indiv.nodes[i]), latent_children = NA, 
                weight = 0)
            names(new_df) = names(coarest_set)
            coarest_set = rbind(coarest_set, new_df)
        }
    }
    return(coarest_set)
}

# Source: functions_v2.R:396
assign_weight=function(tree_df,p,weight.order=-1){
  tree_df$weight[(p+1):nrow(tree_df)]=1
  #print(tree_df)
  nodes=tree_df$node[tree_df$weight!=0]
  result=gather_leaf_nodes_per_non_leaf(tree_df)$df_result
  depth=c(nrow(result))
  for (i in 1:nrow(result)) {
    #depth[i]=1
    #depth[i]=1/sqrt(length(result$leaves[[i]]))
    depth[i]=(length(result$leaves[[i]]))^weight.order
    #depth[i]=1/(length(result$leaves[[i]]))^2
    tree_df$weight[tree_df$node==result$node[i]]=depth[i]
  }
  tree_df$weight=tree_df$weight/mean(depth)
  return(tree_df)
}

# Source: logistic_functions_v2.R:3
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
  #tao=1
  
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
    #iter0=0
    while (!accept){
      #iter0=iter0+1
      eta=Gam-tao*base1/n
      eta.new=prox_tree(eta,lambda = tao*lambda, tree_result)
      if(mean(log(1+exp(crossprod(t(X),eta.new))))-mean(crossprod(t(X),eta.new)*Y) <= mean(log(1+exp.matXGam))-mean(matXGam*Y) + sum((crossprod(X,exp.matXGam/(1+exp.matXGam)-Y)/n)*(eta.new-Gam)) + sum((eta.new-Gam)^2)/(2*tao) || tao <= 1/L0){
        accept=T
      }
      else{
        tao=max(tao/2,1/L0)
      }
    }
    #print(iter0)
    beta0=beta1
    beta1=eta.new
    old.obj=mean(log(1+exp(crossprod(t(X),beta0))))-mean(crossprod(t(X),beta0)*Y)
    new.obj=mean(log(1+exp(crossprod(t(X),beta1))))-mean(crossprod(t(X),beta1)*Y)
    #dis=sqrt(sum((beta1-beta0)^2)/sum(beta0^2))
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

# Source: logistic_functions_v2.R:80
find_max_param.logistic=function (Y, X, coarest_set) {
    p = ncol(X)
    n = nrow(X)
    leaves = 1:p
    #print(length(leaves))
    valid_set = coarest_set[coarest_set$weight != 0, ]
    p1 = p - sum(unlist(lapply(valid_set$leaves, length))) + 
        nrow(valid_set)
    #print(p1)
    Q = matrix(0, nrow = p, ncol = p1)
    remain = rep(1, length(leaves))
    for (i in 1:nrow(valid_set)) {
        ind = match(valid_set[i, ]$leaves[[1]], leaves)
        Q[ind, i] = 1
        remain[ind] = 0
    }
    indiv_num = sum(remain)
    #print(indiv_num)
    if (indiv_num > 0) 
        Q[which(remain == 1), (p1 - indiv_num + 1):p1] = diag(indiv_num)
    X1 = X %*% Q
    fit <- glm(Y ~ X1 + 0, family = binomial(link = "logit"))
    beta1 = coef(fit)
    target = crossprod(X, Y - 1/(1 + exp(-crossprod(t(X1), beta1))))
    vals = sqrt((t(Q) %*% target^2)[1:nrow(valid_set)])/valid_set$weight
    return(max(vals)/n)
}

# Source: logistic_functions_v2.R:136
grid.logistic=function(Y,X,tree_df,true_beta,ridge.param=0,seqc=NULL,thresh=1e-5){
  n=length(Y)
  p=ncol(X)
  tree_result=gather_leaf_nodes_per_non_leaf(tree_df)
  coarest_set=find_coarest(tree_df,tree_result$df_result)
  penalty.max=find_max_param.logistic(Y,X,coarest_set)
  if(is.null(seqc)){
    seqc=exp(seq(-4-log(n),log(penalty.max),length=50))
  }
  #seqc=sort(seqc,decreasing = T)
  #Y.mean=mean(Y)
  Y1=Y
  #X.mean=apply(X, 1, mean)
  X1=X
  #beta0=Y.mean/X.mean
  beta=matrix(0,nrow = p ,ncol = length(seqc))
  #beta[,1]=acc_prox_simple_linear(Y1, X1,tree_result$list_result,lambda = seqc[1],ridge_param = ridge.param,thresh = thresh)$beta1
  beta[,1]=acc_prox_simple_logistic(Y1, X1,tree_result$list_result,lambda = seqc[1],ridge_param = ridge.param,thresh = thresh)
  for (i in 2:length(seqc)) {
    #print(i)
    #beta[,i]=acc_prox_simple_linear(Y1, X1,tree_result$list_result,lambda = seqc[i],warm_start = T,init_beta = as.vector(beta[,i-1]),ridge_param = ridge.param,thresh = thresh)$beta1
    beta[,i]=acc_prox_simple_logistic(Y1, X1,tree_result$list_result,lambda = seqc[i],warm_start = T,init_beta = as.vector(beta[,i-1]),ridge_param = ridge.param,thresh = thresh)
  }
  loss=apply(matrix(rep(true_beta,length(seqc)),nrow = p,byrow = F)-beta,2,function(x) sum(x^2))/p
  return(list(loss=loss,beta=beta,lambda=seqc))
}

# Source: logistic_functions_v2.R:209
negtv_lglkh=function(Y,X,beta){
    n=nrow(X)
    temp=crossprod(t(X),beta)
    return(1/n*as.numeric(-crossprod(Y,temp)+colSums(log(1+exp(temp)))))
}

# Source: logistic_functions_v2.R:335
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

# Source: logistic_functions_v2.R:411
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

# Source: logistic_functions_v2.R:425
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

# Source: logistic_functions_v2.R:458
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

