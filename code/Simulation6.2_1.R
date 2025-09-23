uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

source("/code/functions.R")

nreps=200
n=50
p=100
n0=n*10
ks=seq(10,50,10)
weight.order=c(-1/2,-1,-2)[((uu-1)%/%nreps)+1]

set.seed(uu)
result=data.frame()
for(i in 1:length(ks)){
  single.result=as.data.frame(simulate_error_valid(n=n,p=p,k=ks[i],n0=n0,reps=1,weight.order=weight.order))
  result=rbind(result,single.result)
}

rownames(result)=ks

savename <- paste("/output/Valid_BienTree_n50_p100_kIncre_weightChange/Valid_BienTree_n50_p100_kIncre_weight", weight.order, "_uu", (uu-1)%%nreps + 1, ".RData", sep="")
save(result, file=savename)