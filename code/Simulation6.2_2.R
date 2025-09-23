uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

source("/code/functions.R")

nreps=200
n=500
ps=c(50,100,200,400,600,800,1000)
n0=n*10
rate=0.25
weight.order=c(-1/2,-1,-2)[((uu-1)%/%nreps)+1]

if(sum(((uu-1)%%nreps + 1)==c(74,95,113))){
    set.seed(uu*30)
}else{
    set.seed(uu)
}    
result=data.frame()
for(i in 1:length(ps)){
    p=ps[i]
    k=floor(p*rate)
    k=ifelse(k%%2==0,k,k-1)
    single.result=as.data.frame(simulate_error_valid(n=n,p=p,k=k,n0=n0,reps=1,weight.order=weight.order))
    result=rbind(result,single.result)
    cat(i,',',p,'is finished\n')
}

rownames(result)=ps

savename <- paste("/output/Valid_BienTree_n500_pIncre_kpRatio025_weightChange/Valid_BienTree_n500_pIncre_kpRatio025_weight", weight.order, "_uu", (uu-1)%%nreps + 1, ".RData", sep="")
save(result, file=savename)