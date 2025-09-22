library(Rcpp)
sourceCpp("/code/algos.cpp")
library(dplyr)
library(igraph)
library(MASS)
library(rare)
library(fossil)
library(mclust)
library(profvis)
library(glmnet)


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

################################# This section includes algorithms for squares loss ################################

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

############################ This section simulates the data ################################
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
simulate_data=function(n,group.index,s=0,beta.pre=NULL,ratio=5,new.ratio=NULL,scale.X=FALSE){
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
  #sigma=1
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
simulate_error_tree1_valid=function(n,n0,p0,k,n1=n,s=0,reps=50,ridge.param=0,thresh=1e-5,beta.specify=NULL,tree.specify=NULL,INratio=5,weight.order=-1){
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
    data=simulate_data(n+n0+n1,trees$group,beta.pre = beta.specify,ratio=INratio)
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
