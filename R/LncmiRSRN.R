######################################################################################################
## LncmiRSRN: Identifying lncRNA related miRNA sponge regulatory network
## May 1st, 2017, written by Junpeng Zhang
######################################################################################################

## Extract miRNA-target interactions where targets are from matched lncRNA and mRNA expression data ##
QueryTargetbinding<-function(ExpDataNames,lncR,mR,miRmRbinding,miRlncRbinding){
      
      lncRNames=ExpDataNames[lncR]
      mRNames=ExpDataNames[mR]

      miRmR<-read.csv(miRmRbinding, header=FALSE, sep=",")
      miRmR=as.matrix(miRmR)          

      miRlncR<-read.csv(miRlncRbinding, header=FALSE, sep=",")
      miRlncR=as.matrix(miRlncR)        

      miRmRCandidate=miRmR[which(miRmR[,2] %in% mRNames),]
      miRlncRCandidate=miRlncR[which(miRlncR[,2] %in% lncRNames),]
      return(list(miRlncRCandidate,miRmRCandidate))
}

## Candidate lncRNA-mRNA pairs for miRNA sponge interactions ##
lncRmR<-function(miRlncRCandidate,miRmRCandidate,ExpData,ExpDataNames){
        
        miRlncRCandidate=miRlncRCandidate[sort.list(miRlncRCandidate[,2]),]
        miRlncRCandidate=as.matrix(miRlncRCandidate)
        m1=nrow(miRlncRCandidate)
        n1=ncol(miRlncRCandidate)

        miRmRCandidate=miRmRCandidate[sort.list(miRmRCandidate[,2]),]
        miRmRCandidate=as.matrix(miRmRCandidate)
        m2=nrow(miRmRCandidate)
        n2=ncol(miRmRCandidate)

        lncR=miRlncRCandidate[,2]
        miR1=miRlncRCandidate[,1]
        mR=miRmRCandidate[,2]
        miR2=miRmRCandidate[,1]
        
        lncRNASym=table(lncR)
        mRNASym=table(mR)
        miRSym1=table(miR1)
        miRSym2=table(miR2)

        m3=dim(lncRNASym)
        m4=dim(mRNASym)
        K1=as.matrix(as.numeric(lncRNASym))
        K2=as.matrix(as.numeric(mRNASym))
        lncRNA=as.matrix(rownames(lncRNASym))
        mRNA=as.matrix(rownames(mRNASym))
        
        
        ## Initialize variables
        miRExpIdx=c()
        lncRmRInt=matrix(NA,m3*m4,2)
        R=matrix(NA,m3*m4,5)

        for (i in 1:m3){
              for (j in 1:m4){

                 kk_1=(sum(K1[1:(i-1)])+1):sum(K1[1:i])
                 Interin1=miRlncRCandidate[kk_1,1]
             
                 kk_2=(sum(K2[1:(j-1)])+1):sum(K2[1:j])
                 Interin2=miRmRCandidate[kk_2,1]                 

                 ## Calculate significance of common miRNAs shared by each lncRNA-mRNA pair
		 M1=length(Interin1)
                 M2=length(Interin2)
                 M3=length(intersect(Interin1,Interin2))
                 M4=length(union(rownames(miRSym1),rownames(miRSym2)))
                 M5=1-phyper(M3-1,M2,M4-M2,M1)
                                                      
              if (M3>2 & M5<0.01){
                        
                        lncRmRInt[(i-1)*m3+j,1]=lncRNA[i]
                        lncRmRInt[(i-1)*m3+j,2]=mRNA[j]                      
					
                        lncRExpIdx= which(ExpDataNames==lncRNA[i,1])
                        mRExpIdx= which(ExpDataNames==mRNA[j,1])
                        
			## Calculate Pearson correlation of each lncRNA-mRNA pair
			M6=cor.test(ExpData[,lncRExpIdx],ExpData[,mRExpIdx])$estimate
                        M7=cor.test(ExpData[,lncRExpIdx],ExpData[,mRExpIdx])$p.value
                        
                        R[(i-1)*m3+j,1]=M3
                        R[(i-1)*m3+j,2]=paste(intersect(Interin1,Interin2), collapse = ", ")
		        R[(i-1)*m3+j,3]=M5
                        R[(i-1)*m3+j,4]=M6
			R[(i-1)*m3+j,5]=M7
              }
                                  
          }
        }

## Extract positive correlated lncRNA-mRNA pairs, the p-values are adjusted by BH method.
lncRmRInt=lncRmRInt[which((R[, 1] > 2 & p.adjust(R[, 3],method="BH") < 0.01 & R[, 4] > 0 & p.adjust(R[, 5],method="BH") < 0.01)=='TRUE'),]
RlncRmR=R[which((R[, 1] > 2 & p.adjust(R[, 3],method="BH") < 0.01 & R[, 4] > 0 & p.adjust(R[, 5],method="BH") < 0.01)=='TRUE'),]
PClncRmR=cbind(lncRmRInt, RlncRmR)

return(PClncRmR)
}

## Survival analysis of differential and conserved LncmiRSRN network modules, we use predict.coxph function in survival package to calculate risk scores
SurvAnalyze<-function(ExpData,ExpDataNames,Survival_datacsv,MCLCluster_Name){
library(survival)
library(survcomp)

k=0
hazard_ratio=c()
hazard_ratio_up95=c()
hazard_ratio_low95=c()
myfit=list()
index=c()
p.value=c()
Survival_data=read.csv(Survival_datacsv,header=TRUE,sep=",")
for (i in (1:length(MCLCluster_Name))){
data_Interin=cbind(Survival_data[,2:3],ExpData[,which(ExpDataNames %in% MCLCluster_Name[[i]])])
data_Interin=na.omit(data_Interin)

try_mm=try(coxph(Surv(time, status)~., data = data.frame(data_Interin)),silent=T)
if ('try-error' %in% class(try_mm)) next
mm = coxph(Surv(time, status)~., data = data.frame(data_Interin))

Risk_score=predict(mm,newdata=data.frame(data_Interin),type="risk")

group=rep("NA",dim(data_Interin)[1])
group[Risk_score>median(Risk_score)]="High"
group[Risk_score<=median(Risk_score)]="Low"

data= cbind(data_Interin[,1:2],group)
myfit[[i]]=survfit(Surv(time,status)~group,data=data)

sdf=survdiff(Surv(time,status)~group,data=data)
sdf.p.val=1-pchisq(sdf$chisq, length(sdf$n)-1)
HR = (sdf$obs[1]/sdf$exp[1])/(sdf$obs[2]/sdf$exp[2])
HRup95 = exp(log(HR)+qnorm(0.975)*sqrt(1/sdf$exp[1]+1/sdf$exp[2]))
HRlow95 = exp(log(HR)-qnorm(0.975)*sqrt(1/sdf$exp[1]+1/sdf$exp[2]))

if(sdf.p.val<0.05 & HR>=1.5){
k=k+1
index[k]=i
hazard_ratio[k]=HR
hazard_ratio_up95[k]=HRup95
hazard_ratio_low95[k]=HRlow95
p.value[k]=sdf.p.val
}
}
return(list(index,hazard_ratio,hazard_ratio_up95,hazard_ratio_low95,p.value))
}
