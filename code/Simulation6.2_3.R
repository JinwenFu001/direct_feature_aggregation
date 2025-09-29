uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
source("/code/functions.R")
nreps=200
n=50
ps=c(50,100,200,400,600,800,1000)
n0=n*10
k=20
weight.order=c(-1/2,-1,-2)[((uu-1)%/%nreps)+1]

set.seed(uu)
result=data.frame()
for(i in 1:length(ps)){
    p=ps[i]
    single.result=as.data.frame(simulate_error_valid.logistic(n=n,p=p,k=k,n0=n0,reps=1,weight.order=weight.order))
    result=rbind(result,single.result)
    cat(i,',',p,'is finished\n')
}

rownames(result)=ps

savename <- paste("/output/Valid_BienTree_Logistic_n50_pIncre_k20_weightChange/Valid_BienTree_Logistic_n50_pIncre_k20_weight", weight.order, "_uu", (uu-1)%%nreps + 1, ".RData", sep="")
save(result, file=savename)