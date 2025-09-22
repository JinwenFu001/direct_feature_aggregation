uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

source("/code/functions.R")

nreps=200
n=50
p0s=c(3,5,10,20,30,40,50)
n0=n*10
k=20
weight.order=c(-1/2,-1,-2)[((uu-1)%/%nreps)+1]

if(sum(uu==c(55,131,126,146,153))) uu=uu*30
set.seed(uu)
result=data.frame()
for(i in 1:length(p0s)){
    p0=p0s[i]
    single.result=as.data.frame(simulate_error_tree1_valid(n=n,n0=n0,p0=p0,k=k,reps=1,weight.order=weight.order))
    result=rbind(result,single.result)
    cat(i,',',p0,'is finished\n')
}

rownames(result)=p0s

savename <- paste("/output/Valid_Tree1_n50_p0Incre_k20_scaled_weightChange/Valid_Tree1_n50_p0Incre_k20_scaled_weight", weight.order, "_uu", (uu-1)%%nreps + 1, ".RData", sep="")
save(result, file=savename)