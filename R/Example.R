#############################################################################
## Main scripts for identifying lncRNA related miRNA sponge regulatory network ## 
#############################################################################
source("LncmiRSRN.R")

## Import matched lncRNA and mRNA expression data
load("Expression_Data.RData")

## Extracting candidate miRNA-target interactions by integrating putative miRNA-target interactions with expression data
miRmRbinding="miRTarBase_v7.0+TarBase_v7.0_High.csv"
miRlncRbinding="NPInter_v3.0+LncBase_v2.csv"
lncR=1:9704
mR=9705:27986
lncRNames=ExpDataNames[lncR]
mRNames=ExpDataNames[mR]
miRTargetCandidate=QueryTargetbinding(ExpDataNames,lncR,mR,miRmRbinding,miRlncRbinding)
miRlncRCandidate=miRTargetCandidate[[1]]
miRmRCandidate=miRTargetCandidate[[2]]

## Identifying lncRNA-mRNA related miRNA sponge interactions in GBM, LSCC, OvCa and PrCa
PClncRmR_GBM=lncRmR(miRlncRCandidate,miRmRCandidate,ExpData_GBM,ExpDataNames)
PClncRmR_LSCC=lncRmR(miRlncRCandidate,miRmRCandidate,ExpData_LSCC,ExpDataNames)
PClncRmR_OvCa=lncRmR(miRlncRCandidate,miRmRCandidate,ExpData_OvCa,ExpDataNames)
PClncRmR_PrCa=lncRmR(miRlncRCandidate,miRmRCandidate,ExpData_PrCa,ExpDataNames)

library("igraph")
PClncRmR_GBM_graph=make_graph(c(t(PClncRmR_GBM[,1:2])),directed = TRUE)
PClncRmR_LSCC_graph=make_graph(c(t(PClncRmR_LSCC[,1:2])),directed = TRUE)
PClncRmR_OvCa_graph=make_graph(c(t(PClncRmR_OvCa[,1:2])),directed = TRUE)
PClncRmR_PrCa_graph=make_graph(c(t(PClncRmR_PrCa[,1:2])),directed = TRUE)

## By using ParallelPC R package, we generate lncRNA-mRNA regulatory relationships from matched lncRNA and mRNA expression data. The lncRNAs and mRNAs are from lncRNA-mRNA related miRNA sponge interactions
library(ParallelPC)
library(bnlearn)
library(pcalg)
library(parallel)
Sponge_Exp_GBM=ExpData_GBM[,c(which(ExpDataNames %in% unique(PClncRmR_GBM[,1])),which(ExpDataNames %in% unique(PClncRmR_GBM[,2])))]
Sponge_Exp_LSCC=ExpData_LSCC[,c(which(ExpDataNames %in% unique(PClncRmR_LSCC[,1])),which(ExpDataNames %in% unique(PClncRmR_LSCC[,2])))]
Sponge_Exp_OvCa=ExpData_OvCa[,c(which(ExpDataNames %in% unique(PClncRmR_OvCa[,1])),which(ExpDataNames %in% unique(PClncRmR_OvCa[,2])))]
Sponge_Exp_PrCa=ExpData_PrCa[,c(which(ExpDataNames %in% unique(PClncRmR_PrCa[,1])),which(ExpDataNames %in% unique(PClncRmR_PrCa[,2])))]

# Using ParallelPC R package to estimate causal effects between sponge lncRNAs and mRNAs
Causal_Score_GBM=IDA_parallel(Sponge_Exp_GBM,1:length(unique(PClncRmR_GBM[,1])),(length(unique(PClncRmR_GBM[,1]))+1):length(c(unique(PClncRmR_GBM[,1]),unique(PClncRmR_GBM[,2]))),"parallel",0.01, 6)
Causal_Score_LSCC=IDA_parallel(Sponge_Exp_LSCC,1:length(unique(PClncRmR_LSCC[,1])),(length(unique(PClncRmR_LSCC[,1]))+1):length(c(unique(PClncRmR_LSCC[,1]),unique(PClncRmR_LSCC[,2]))),"parallel",0.01, 6)
Causal_Score_OvCa=IDA_parallel(Sponge_Exp_OvCa,1:length(unique(PClncRmR_OvCa[,1])),(length(unique(PClncRmR_OvCa[,1]))+1):length(c(unique(PClncRmR_OvCa[,1]),unique(PClncRmR_OvCa[,2]))),"parallel",0.01, 6)
Causal_Score_PrCa=IDA_parallel(Sponge_Exp_PrCa,1:length(unique(PClncRmR_PrCa[,1])),(length(unique(PClncRmR_PrCa[,1]))+1):length(c(unique(PClncRmR_PrCa[,1]),unique(PClncRmR_PrCa[,2]))),"parallel",0.01, 6)

# Calculate the Fisher's asymptotic p-values based on the causal effects that sponge lncRNAs on mRNAs 
library(WGCNA)
Causal_Score_GBM_normalize=Causal_Score_GBM
Causal_Score_LSCC_normalize=Causal_Score_LSCC
Causal_Score_OvCa_normalize=Causal_Score_OvCa
Causal_Score_PrCa_normalize=Causal_Score_PrCa

Causal_Score_GBM_normalize[which(Causal_Score_GBM_normalize>=1,arr.ind = TRUE)]=1
Causal_Score_GBM_normalize[which(Causal_Score_GBM_normalize<=-1,arr.ind = TRUE)]=-1
Causal_Score_LSCC_normalize[which(Causal_Score_LSCC_normalize>=1,arr.ind = TRUE)]=1
Causal_Score_LSCC_normalize[which(Causal_Score_LSCC_normalize<=-1,arr.ind = TRUE)]=-1
Causal_Score_OvCa_normalize[which(Causal_Score_OvCa_normalize>=1,arr.ind = TRUE)]=1
Causal_Score_OvCa_normalize[which(Causal_Score_OvCa_normalize<=-1,arr.ind = TRUE)]=-1
Causal_Score_PrCa_normalize[which(Causal_Score_PrCa_normalize>=1,arr.ind = TRUE)]=1
Causal_Score_PrCa_normalize[which(Causal_Score_PrCa_normalize<=-1,arr.ind = TRUE)]=-1

Causal_Score_GBM_pvalue=corPvalueFisher(Causal_Score_GBM_normalize, dim(Sponge_Exp_GBM)[1])
Causal_Score_LSCC_pvalue=corPvalueFisher(Causal_Score_LSCC_normalize, dim(Sponge_Exp_LSCC)[1])
Causal_Score_OvCa_pvalue=corPvalueFisher(Causal_Score_OvCa_normalize, dim(Sponge_Exp_OvCa)[1])
Causal_Score_PrCa_pvalue=corPvalueFisher(Causal_Score_PrCa_normalize, dim(Sponge_Exp_PrCa)[1])

Causal_Score_GBM_pvalue_adjust=matrix(p.adjust(c(Causal_Score_GBM_pvalue), method="BH"), nrow=dim(Causal_Score_GBM_pvalue)[1])
Causal_Score_LSCC_pvalue_adjust=matrix(p.adjust(c(Causal_Score_LSCC_pvalue), method="BH"), nrow=dim(Causal_Score_LSCC_pvalue)[1])
Causal_Score_OvCa_pvalue_adjust=matrix(p.adjust(c(Causal_Score_OvCa_pvalue), method="BH"), nrow=dim(Causal_Score_OvCa_pvalue)[1])
Causal_Score_PrCa_pvalue_adjust=matrix(p.adjust(c(Causal_Score_PrCa_pvalue), method="BH"), nrow=dim(Causal_Score_PrCa_pvalue)[1])

## Using Fisher's asymptotic p-values to evaluate strength of causal effects, the p-value cutoff is 0.05
Causal_lncRmR_GBM=cbind(colnames(Causal_Score_GBM_pvalue)[which(Causal_Score_GBM_pvalue_adjust<0.05,arr.ind = TRUE)[,2]],rownames(Causal_Score_GBM_pvalue)[which(Causal_Score_GBM_pvalue_adjust<0.05,arr.ind = TRUE)[,1]])
Causal_lncRmR_LSCC=cbind(colnames(Causal_Score_LSCC_pvalue)[which(Causal_Score_LSCC_pvalue_adjust<0.05,arr.ind = TRUE)[,2]],rownames(Causal_Score_LSCC_pvalue)[which(Causal_Score_LSCC_pvalue_adjust<0.05,arr.ind = TRUE)[,1]])
Causal_lncRmR_OvCa=cbind(colnames(Causal_Score_OvCa_pvalue)[which(Causal_Score_OvCa_pvalue_adjust<0.05,arr.ind = TRUE)[,2]],rownames(Causal_Score_OvCa_pvalue)[which(Causal_Score_OvCa_pvalue_adjust<0.05,arr.ind = TRUE)[,1]])
Causal_lncRmR_PrCa=cbind(colnames(Causal_Score_PrCa_pvalue)[which(Causal_Score_PrCa_pvalue_adjust<0.05,arr.ind = TRUE)[,2]],rownames(Causal_Score_PrCa_pvalue)[which(Causal_Score_PrCa_pvalue_adjust<0.05,arr.ind = TRUE)[,1]])

Causal_Score_GBM_Effect=Causal_Score_GBM[Causal_Score_GBM_pvalue_adjust<0.05]
Causal_Score_LSCC_Effect=Causal_Score_LSCC[Causal_Score_LSCC_pvalue_adjust<0.05]
Causal_Score_OvCa_Effect=Causal_Score_OvCa[Causal_Score_OvCa_pvalue_adjust<0.05]
Causal_Score_PrCa_Effect=Causal_Score_PrCa[Causal_Score_PrCa_pvalue_adjust<0.05]

## Generate lncRNA-related miRNA sponge regulatory networks
library("igraph")
Causal_lncRmR_GBM_graph=make_graph(c(t(Causal_lncRmR_GBM)),directed = TRUE)
Causal_lncRmR_LSCC_graph=make_graph(c(t(Causal_lncRmR_LSCC)),directed = TRUE)
Causal_lncRmR_OvCa_graph=make_graph(c(t(Causal_lncRmR_OvCa)),directed = TRUE)
Causal_lncRmR_PrCa_graph=make_graph(c(t(Causal_lncRmR_PrCa)),directed = TRUE)

Causal_Sponge_lncRmR_GBM_graph=Causal_lncRmR_GBM_graph %s% PClncRmR_GBM_graph
Causal_Sponge_lncRmR_LSCC_graph=Causal_lncRmR_LSCC_graph %s% PClncRmR_LSCC_graph
Causal_Sponge_lncRmR_OvCa_graph=Causal_lncRmR_OvCa_graph %s% PClncRmR_OvCa_graph
Causal_Sponge_lncRmR_PrCa_graph=Causal_lncRmR_PrCa_graph %s% PClncRmR_PrCa_graph

# Degree distribution of lncRNA-related miRNA sponge regulatory networks
Causal_Sponge_lncRmR_GBM_Degree=degree(graph_from_data_frame(as_data_frame(Causal_Sponge_lncRmR_GBM_graph),directed=FALSE))
Causal_Sponge_lncRmR_LSCC_Degree=degree(graph_from_data_frame(as_data_frame(Causal_Sponge_lncRmR_LSCC_graph),directed=FALSE))
Causal_Sponge_lncRmR_OvCa_Degree=degree(graph_from_data_frame(as_data_frame(Causal_Sponge_lncRmR_OvCa_graph),directed=FALSE))
Causal_Sponge_lncRmR_PrCa_Degree=degree(graph_from_data_frame(as_data_frame(Causal_Sponge_lncRmR_PrCa_graph),directed=FALSE))

## hub lncRNAs in GBM, LSCC, OvCa and PrCa dataset
Causal_Sponge_lncR_GBM_outdegree=degree(Causal_Sponge_lncRmR_GBM_graph,mode="out")
hub_lncRNAs_GBM=names(sort(Causal_Sponge_lncR_GBM_outdegree[which(Causal_Sponge_lncR_GBM_outdegree!=0)],decreasing=TRUE))[1:ceiling(0.2*length(which(Causal_Sponge_lncR_GBM_outdegree!=0)))]
Causal_Sponge_lncR_LSCC_outdegree=degree(Causal_Sponge_lncRmR_LSCC_graph,mode="out")
hub_lncRNAs_LSCC=names(sort(Causal_Sponge_lncR_LSCC_outdegree[which(Causal_Sponge_lncR_LSCC_outdegree!=0)],decreasing=TRUE))[1:ceiling(0.2*length(which(Causal_Sponge_lncR_LSCC_outdegree!=0)))]
Causal_Sponge_lncR_OvCa_outdegree=degree(Causal_Sponge_lncRmR_OvCa_graph,mode="out")
hub_lncRNAs_OvCa=names(sort(Causal_Sponge_lncR_OvCa_outdegree[which(Causal_Sponge_lncR_OvCa_outdegree!=0)],decreasing=TRUE))[1:ceiling(0.2*length(which(Causal_Sponge_lncR_OvCa_outdegree!=0)))]
Causal_Sponge_lncR_PrCa_outdegree=degree(Causal_Sponge_lncRmR_PrCa_graph,mode="out")
hub_lncRNAs_PrCa=names(sort(Causal_Sponge_lncR_PrCa_outdegree[which(Causal_Sponge_lncR_PrCa_outdegree!=0)],decreasing=TRUE))[1:ceiling(0.2*length(which(Causal_Sponge_lncR_PrCa_outdegree!=0)))]

## Similarity between GBM, LSCC, OvCa and PrCa in terms of hub lncRNAs
GBM_LSCC_hub_Similarity=length(intersect(hub_lncRNAs_GBM,hub_lncRNAs_LSCC))/min(length(hub_lncRNAs_GBM),length(hub_lncRNAs_LSCC))
GBM_OvCa_hub_Similarity=length(intersect(hub_lncRNAs_GBM,hub_lncRNAs_OvCa))/min(length(hub_lncRNAs_GBM),length(hub_lncRNAs_OvCa))
GBM_PrCa_hub_Similarity=length(intersect(hub_lncRNAs_GBM,hub_lncRNAs_PrCa))/min(length(hub_lncRNAs_GBM),length(hub_lncRNAs_PrCa))
LSCC_OvCa_hub_Similarity=length(intersect(hub_lncRNAs_LSCC,hub_lncRNAs_OvCa))/min(length(hub_lncRNAs_LSCC),length(hub_lncRNAs_OvCa))
LSCC_PrCa_hub_Similarity=length(intersect(hub_lncRNAs_LSCC,hub_lncRNAs_PrCa))/min(length(hub_lncRNAs_LSCC),length(hub_lncRNAs_PrCa))
OvCa_PrCa_hub_Similarity=length(intersect(hub_lncRNAs_OvCa,hub_lncRNAs_PrCa))/min(length(hub_lncRNAs_OvCa),length(hub_lncRNAs_PrCa))

## Differential and conserved LncmiRSRNs between GBM, LSCC, OvCa and PrCa
Causal_Sponge_lncRmR_Unique_graph=(Causal_Sponge_lncRmR_GBM_graph %u% Causal_Sponge_lncRmR_LSCC_graph) %u% (Causal_Sponge_lncRmR_OvCa_graph %u% Causal_Sponge_lncRmR_PrCa_graph)
Causal_Sponge_lncRmR_Unique_Interactions=as_data_frame(Causal_Sponge_lncRmR_Unique_graph)

Causal_Sponge_lncRmR_Differential_graph=(Causal_Sponge_lncRmR_GBM_graph %m% Causal_Sponge_lncRmR_LSCC_graph %m% Causal_Sponge_lncRmR_OvCa_graph %m% Causal_Sponge_lncRmR_PrCa_graph) %u% 
                                        (Causal_Sponge_lncRmR_LSCC_graph %m% Causal_Sponge_lncRmR_GBM_graph %m% Causal_Sponge_lncRmR_OvCa_graph %m% Causal_Sponge_lncRmR_PrCa_graph) %u% 
                                        (Causal_Sponge_lncRmR_OvCa_graph %m% Causal_Sponge_lncRmR_GBM_graph %m% Causal_Sponge_lncRmR_LSCC_graph %m% Causal_Sponge_lncRmR_PrCa_graph) %u% 
                                        (Causal_Sponge_lncRmR_PrCa_graph %m% Causal_Sponge_lncRmR_GBM_graph %m% Causal_Sponge_lncRmR_LSCC_graph %m% Causal_Sponge_lncRmR_OvCa_graph)
Causal_Sponge_lncRmR_Differential_Interactions=as_data_frame(Causal_Sponge_lncRmR_Differential_graph)

Causal_Sponge_lncRmR_Conserved_graph=Causal_Sponge_lncRmR_Unique_graph %m% Causal_Sponge_lncRmR_Differential_graph
Causal_Sponge_lncRmR_Conserved_Interactions=as_data_frame(Causal_Sponge_lncRmR_Conserved_graph)

# Degree distribution of differential and conserved LncmiRSRNs
Causal_Sponge_lncRmR_Differential_Degree=degree(graph_from_data_frame(as_data_frame(Causal_Sponge_lncRmR_Differential_graph),directed=FALSE))
Causal_Sponge_lncRmR_Conserved_Degree=degree(graph_from_data_frame(as_data_frame(Causal_Sponge_lncRmR_Conserved_graph),directed=FALSE))

## Similarity between GBM, LSCC, OvCa and PrCa in terms of sponge lncRNA-mRNA regulatory relationships
GBM_LSCC_Network_Similarity=dim(as_data_frame(Causal_Sponge_lncRmR_GBM_graph %s% Causal_Sponge_lncRmR_LSCC_graph))[1]/min(dim(as_data_frame(Causal_Sponge_lncRmR_GBM_graph))[1],dim(as_data_frame(Causal_Sponge_lncRmR_LSCC_graph))[1])
GBM_OvCa_Network_Similarity=dim(as_data_frame(Causal_Sponge_lncRmR_GBM_graph %s% Causal_Sponge_lncRmR_OvCa_graph))[1]/min(dim(as_data_frame(Causal_Sponge_lncRmR_GBM_graph))[1],dim(as_data_frame(Causal_Sponge_lncRmR_OvCa_graph))[1])
GBM_PrCa_Network_Similarity=dim(as_data_frame(Causal_Sponge_lncRmR_GBM_graph %s% Causal_Sponge_lncRmR_PrCa_graph))[1]/min(dim(as_data_frame(Causal_Sponge_lncRmR_GBM_graph))[1],dim(as_data_frame(Causal_Sponge_lncRmR_PrCa_graph))[1])
LSCC_OvCa_Network_Similarity=dim(as_data_frame(Causal_Sponge_lncRmR_LSCC_graph %s% Causal_Sponge_lncRmR_OvCa_graph))[1]/min(dim(as_data_frame(Causal_Sponge_lncRmR_LSCC_graph))[1],dim(as_data_frame(Causal_Sponge_lncRmR_OvCa_graph))[1])
LSCC_PrCa_Network_Similarity=dim(as_data_frame(Causal_Sponge_lncRmR_LSCC_graph %s% Causal_Sponge_lncRmR_PrCa_graph))[1]/min(dim(as_data_frame(Causal_Sponge_lncRmR_LSCC_graph))[1],dim(as_data_frame(Causal_Sponge_lncRmR_PrCa_graph))[1])
OvCa_PrCa_Network_Similarity=dim(as_data_frame(Causal_Sponge_lncRmR_OvCa_graph %s% Causal_Sponge_lncRmR_PrCa_graph))[1]/min(dim(as_data_frame(Causal_Sponge_lncRmR_OvCa_graph))[1],dim(as_data_frame(Causal_Sponge_lncRmR_PrCa_graph))[1])

## Identify differential and conserved LncmiRSRN network modules using MCL method
library("ProNet")
Causal_Sponge_lncRmR_Differential_MCLCluster=cluster(graph_from_data_frame(Causal_Sponge_lncRmR_Differential_Interactions,directed=FALSE),method="MCL",expansion = 2,inflation = 2,directed = FALSE,layout="fruchterman.reingold")
res=Causal_Sponge_lncRmR_Differential_MCLCluster
fileName="Causal_Sponge_lncRmR_Differential_MCLCluster.txt"
MCLCluster_Differential_Name=list()
k1=0
for (i in 1:max(res)) {
    lncR=length(which(rownames(as.matrix(res))[which(res==i)] %in% lncRNames))
    mR=length(which(rownames(as.matrix(res))[which(res==i)] %in% mRNames))
  if(length(which(res==i))>=4 & lncR>1 & mR>1)
  {
    k1=k1+1
    MCLCluster_Differential_Name[[k1]]=rownames(as.matrix(res))[which(res==i)]    
    cat(c(k1,'\t',length(which(res==i))),'\t',lncR,'\t',mR,file=fileName,sep="",append=TRUE)
    for(j in which(res==i))
    {
      cat(c('\t',rownames(as.matrix(res))[j]),file=fileName,sep="",append=TRUE)
    }
    if(i!=max(res))
    {
      cat('\n',file=fileName,sep="",append=TRUE)
    }
  }
} 

Causal_Sponge_lncRmR_Conserved_MCLCluster=cluster(graph_from_data_frame(Causal_Sponge_lncRmR_Conserved_Interactions,directed=FALSE),method="MCL",expansion = 2,inflation = 2,directed = FALSE,layout="fruchterman.reingold")
res=Causal_Sponge_lncRmR_Conserved_MCLCluster
fileName="Causal_Sponge_lncRmR_Conserved_MCLCluster.txt"
MCLCluster_Conserved_Name=list()
k2=0
for (i in 1:max(res)) {
    lncR=length(which(rownames(as.matrix(res))[which(res==i)] %in% lncRNames))
    mR=length(which(rownames(as.matrix(res))[which(res==i)] %in% mRNames))
  if(length(which(res==i))>=4 & lncR>1 & mR>1)
  {
    k2=k2+1
    MCLCluster_Conserved_Name[[k2]]=rownames(as.matrix(res))[which(res==i)]
    cat(c(k2,'\t',length(which(res==i))),'\t',lncR,'\t',mR,file=fileName,sep="",append=TRUE)
    for(j in which(res==i))
    {
      cat(c('\t',rownames(as.matrix(res))[j]),file=fileName,sep="",append=TRUE)
    }
    if(i!=max(res))
    {
      cat('\n',file=fileName,sep="",append=TRUE)
    }
  }
}
 

## Survival analysis of differential and conserved LncmiRSRN network modules, we use predict.coxph function in survival package to calculate risk scores
GBM_Differential_Survival=SurvAnalyze(ExpData_GBM,ExpDataNames,"GBM_survival.csv",MCLCluster_Differential_Name)
LSCC_Differential_Survival=SurvAnalyze(ExpData_LSCC,ExpDataNames,"LSCC_survival.csv",MCLCluster_Differential_Name)
OvCa_Differential_Survival=SurvAnalyze(ExpData_OvCa,ExpDataNames,"OvCa_survival.csv",MCLCluster_Differential_Name)
PrCa_Differential_Survival=SurvAnalyze(ExpData_PrCa,ExpDataNames,"PrCa_survival.csv",MCLCluster_Differential_Name)

GBM_Conserved_Survival=SurvAnalyze(ExpData_GBM,ExpDataNames,"GBM_survival.csv",MCLCluster_Conserved_Name)
LSCC_Conserved_Survival=SurvAnalyze(ExpData_LSCC,ExpDataNames,"LSCC_survival.csv",MCLCluster_Conserved_Name)
OvCa_Conserved_Survival=SurvAnalyze(ExpData_OvCa,ExpDataNames,"OvCa_survival.csv",MCLCluster_Conserved_Name)
PrCa_Conserved_Survival=SurvAnalyze(ExpData_PrCa,ExpDataNames,"PrCa_survival.csv",MCLCluster_Conserved_Name)

## GO and KEGG enrichment analysis of differential and conserved LncmiRSRN network modules
library(clusterProfiler)
library(org.Hs.eg.db)

Differential_entrezIDs=lapply(1:length(MCLCluster_Differential_Name), function(i) mget(MCLCluster_Differential_Name[[i]], org.Hs.egSYMBOL2EG, ifnotfound=NA))
Differential_entrezIDs=lapply(1:length(MCLCluster_Differential_Name), function(i) as.character(Differential_entrezIDs[[i]]))
Differential_enrichGO=lapply(1:length(MCLCluster_Differential_Name), function(i) enrichGO(Differential_entrezIDs[[i]],OrgDb='org.Hs.eg.db',ont = "BP",pvalueCutoff=0.05))
Differential_enrichKEGG=lapply(1:length(MCLCluster_Differential_Name), function(i) enrichKEGG(Differential_entrezIDs[[i]], organism="hsa",pvalueCutoff=0.05))

Conserved_entrezIDs=lapply(1:length(MCLCluster_Conserved_Name), function(i) mget(MCLCluster_Conserved_Name[[i]], org.Hs.egSYMBOL2EG, ifnotfound=NA))
Conserved_entrezIDs=lapply(1:length(MCLCluster_Conserved_Name), function(i) as.character(Conserved_entrezIDs[[i]]))
Conserved_enrichGO=lapply(1:length(MCLCluster_Conserved_Name), function(i) enrichGO(Conserved_entrezIDs[[i]],OrgDb='org.Hs.eg.db',ont = "BP",pvalueCutoff=0.05))
Conserved_enrichKEGG=lapply(1:length(MCLCluster_Conserved_Name), function(i) enrichKEGG(Conserved_entrezIDs[[i]], organism="hsa",pvalueCutoff=0.05))

## Cancer gene enrichment of differential and conserved LncmiRSRN network modules
GBM_genes=read.csv("GBM_genes.csv",header=FALSE)
LSCC_genes=read.csv("LSCC_genes.csv",header=FALSE)
OvCa_genes=read.csv("OvCa_genes.csv",header=FALSE)
PrCa_genes=read.csv("PrCa_genes.csv",header=FALSE)

Differential_GBM_Num=lapply(1:length(MCLCluster_Differential_Name), function(i) length(intersect(MCLCluster_Differential_Name[[i]],c(as.matrix(GBM_genes)))))
Differential_LSCC_Num=lapply(1:length(MCLCluster_Differential_Name), function(i) length(intersect(MCLCluster_Differential_Name[[i]],c(as.matrix(LSCC_genes)))))
Differential_OvCa_Num=lapply(1:length(MCLCluster_Differential_Name), function(i) length(intersect(MCLCluster_Differential_Name[[i]],c(as.matrix(OvCa_genes)))))
Differential_PrCa_Num=lapply(1:length(MCLCluster_Differential_Name), function(i) length(intersect(MCLCluster_Differential_Name[[i]],c(as.matrix(PrCa_genes)))))

Conserved_GBM_Num=lapply(1:length(MCLCluster_Conserved_Name), function(i) length(intersect(MCLCluster_Conserved_Name[[i]],c(as.matrix(GBM_genes)))))
Conserved_LSCC_Num=lapply(1:length(MCLCluster_Conserved_Name), function(i) length(intersect(MCLCluster_Conserved_Name[[i]],c(as.matrix(LSCC_genes)))))
Conserved_OvCa_Num=lapply(1:length(MCLCluster_Conserved_Name), function(i) length(intersect(MCLCluster_Conserved_Name[[i]],c(as.matrix(OvCa_genes)))))
Conserved_PrCa_Num=lapply(1:length(MCLCluster_Conserved_Name), function(i) length(intersect(MCLCluster_Conserved_Name[[i]],c(as.matrix(PrCa_genes)))))

save.image("LncmiRSRN_p.adjusted_value=0.01_GBM+LSCC+OvCa+PrCa.RData")
