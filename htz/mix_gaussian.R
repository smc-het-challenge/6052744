options(show.error.messages=F,warn=0)
args = commandArgs(trailingOnly = T)
raw_vcf = read.table(args[1],sep="\t",comment.char = "#")
data_col = 11
name_col = 9
attribute_names = unlist(strsplit(as.vector(raw_vcf[1,name_col]),":"))
AD_ind = which(attribute_names=="AD")
tumor_data = lapply(as.vector(raw_vcf[,data_col]), function(x){ unlist(strsplit(x,":")) })
tumor_stats_AD = lapply(tumor_data,function(x){ unlist(strsplit(x[AD_ind],",")) })
mutation_reads = unlist(lapply(tumor_stats_AD,function(x){as.numeric(x[2])}))
non_malignant_reads = unlist(lapply(tumor_stats_AD,function(x){as.numeric(x[1])}))
total_reads = mutation_reads+non_malignant_reads
mutation_rate = mutation_reads/total_reads

a=suppressMessages(require(mixtools))
if(!a)
{
  install.packages("mixtools")
  suppressMessages(require(mixtools))
}


tmp = kmeans(mutation_rate,centers = 2)
ss=tmp$tot.withinss/tmp$totss
{if(ss<0.1)
{
  mu1 = mean(mutation_rate)
  sd1 = sd(mutation_rate)
  cellarity = mu1
  No_cluster = 1
  counts = length(mutation_rate)
  result_1c = c(1,counts,mu1)
  assignment_2a = rep(1,counts)
  ccm_2b = matrix(1, ncol=counts,nrow=counts)
  tree_3a = c(1,0)
  adm_3b=matrix(0,ncol=counts,nrow=counts)
}
else
{
  tmp = lapply(1:9,function(x){suppressMessages(normalmixEM(mutation_rate,k= (x+1)))})
  loglikhood = unlist(lapply(tmp,function(x){x$loglik}))
  penal = (2*(2:10)+1:9)*log(length(mutation_rate))
  bics = -loglikhood+penal
  bic_diff = (bics[2:9]-bics[1:8])/bics[1:8]
  {if(length(which(bic_diff[1:min(which(bic_diff>0))]>0.01))>0)
  {
    No_cluster = max(which(bic_diff[1:min(which(bic_diff>0))] > 0.01))+2
  }
  else
  {
    No_cluster = 2
  }}
  mix_results = tmp[[No_cluster-1]]
  cellarity = max(mix_results$mu)
  assignment_2a = apply(mix_results$posterior,1,which.max)
  counts = unlist(lapply(1:No_cluster,function(x){length(which(assignment_2a==x))}))
  result_1c = cbind(1:No_cluster,counts,mix_results$mu)
  sort.mu = sort(mix_results$mu,decreasing = T)
  tree_3a = NULL
  for(i in 1:No_cluster)
  {
    tree_3a=rbind(tree_3a,c(i,which(sort.mu==mix_results$mu[i])-1))
  }
  ccm_2b = diag(1,nrow = length(mutation_rate))
  adm_3b=matrix(0,nrow = length(mutation_rate),ncol = length(mutation_rate))
  for(i in 1:(length(mutation_rate)-1))
  {
    for(j in (i+1):length(mutation_rate))
    {
     ccm_2b[i,j] = ccm_2b[j,i] = round(sum(mix_results$posterior[i,]*mix_results$posterior[j,]),4)
     if(assignment_2a[i]==assignment_2a[j])
     {
       adm_3b[i,j]=adm_3b[j,i]=0
     }
     else
     {
       adm_3b[i,j] = as.numeric(tree_3a[assignment_2a[i],2]<tree_3a[assignment_2a[j],2])*(1-ccm_2b[i,j])
       adm_3b[j,i] = as.numeric(tree_3a[assignment_2a[j],2]<tree_3a[assignment_2a[i],2])*(1-ccm_2b[j,i])
     }
    }
    
  }
  
}}

write.table(cellarity,"subchallenge1A.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(No_cluster,"subchallenge1B.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(result_1c,"subchallenge1C.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(assignment_2a,"subchallenge2A.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(ccm_2b,"subchallenge2B.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(tree_3a,"subchallenge3A.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(adm_3b,"subchallenge3B.txt",row.names=F,col.names=F,quote=F,sep="\t")
