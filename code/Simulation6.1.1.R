uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

source("/code/scaled_functions.R")
source("/code/vary_tree.R")
set.seed(20)
n=50
n0=500
k=10
p=60
nreps=200
weight.order=c(-1/2,-1,-2)[((uu-1)%/%nreps)+1]
trees=simulate_tree(k = k,p = p)
tree_df0=hclust_to_df(trees$tree,weight.order)
merge.nodes=unlist(find_latent_nodes_with_exact_leaf_groups(tree_df0 ,trees$group))
L1=find_all_parents(tree_df0,merge.nodes)
L0=(1:nrow(tree_df0))[-c(1:p,L1,merge.nodes)]
l0_seq=c(10,20,30)
l1_seq=c(3,6)
tree1_l0=tree_df0; tree1_l0$weight[sample(L0,l0_seq[1],replace = F)]=0
tree2_l0=tree_df0; tree2_l0$weight[sample(L0,l0_seq[2],replace = F)]=0
tree3_l0=tree_df0; tree3_l0$weight[sample(L0,l0_seq[3],replace = F)]=0
tree1_l1=tree_df0; tree1_l1$weight[sample(L1,l1_seq[1],replace = F)]=0
tree2_l1=tree_df0; tree2_l1$weight[sample(L1,l1_seq[2],replace = F)]=0

df_list=list(tree_df0,tree1_l0,tree2_l0,tree3_l0,tree1_l1,tree2_l1)



result=data.frame()
for(i in 1:length(df_list)){
    set.seed(uu)
    tree_df=df_list[[i]]
    single.result=as.data.frame(simulate_error_given_tree(n=n,n0=n0, tree_df=tree_df,tree_group=trees$group,n1=n,s=0,reps=1,ridge.param=0,thresh=1e-5,INratio=5))
    result=rbind(result,single.result)
    cat(i,'is finished\n')
}

rownames(result)=c("tree_df0","tree1_l0","tree2_l0","tree3_l0","tree1_l1","tree2_l1")

savename <- paste("output/Valid_BienTree_n50_p60_k10_treeChange_weightChange/Valid_BienTree_n50_p60_k10_treeChange_weight", weight.order, "_uu", (uu-1)%%nreps + 1, ".RData", sep="")
save(result, file=savename)