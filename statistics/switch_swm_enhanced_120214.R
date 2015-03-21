library(linkcomm)
library(WGCNA)
library(psych)
library(R.matlab)
library(e1071)
library(calibrate)
library(gplot)
library(sciplot)


motion<- as.matrix(read.csv(file="/scr/melisse1/NKI_enhanced/nki-rs_lite_r1_v1_RfMRI_mx_645_rest_motion_parameters.csv",head=TRUE,sep=",")) #motion information
pheno <- as.matrix(read.csv(file="/scr/melisse1/NKI_enhanced/nki-rs_lite_r1_phenotypic_v1.csv",head=TRUE,sep=",")) # this file comes from NKI but is not open to the public

### find subjects with motion parameters
pheno_new=matrix(nrow=dim(motion)[1],ncol=dim(pheno)[2])
subject_ids=as.integer(rownames(motion))
count=1;
for (i in (1:length(subject_ids))){
    idx=which(pheno[,1]==subject_ids[i]);
    if(length(idx)>0){
        print(count)
        pheno_new[count,]=as.matrix(pheno[idx,]);
        count=count+1;
    }
}
pheno=pheno_new
## exclude corrupted images
idx_corr=which(pheno[,1]=='165660' | pheno[,1]=='181439') # remvoe subjects with imaging problems
pheno=pheno[-idx_corr,]
motion=motion[-idx_corr,]

motion_max=motion[,8]
## exclude first time
idx_max=which(motion_max>3) # remove subjects over max motion
pheno=pheno[-idx_max,]
motion=motion[-idx_max,]

as.numeric(pheno[1:4,2])
motion_max=as.numeric(motion[,3])
motion_fd=as.numeric(motion[,33])

## exclude second time
max_fd=(mean(motion_fd))
idx_fd=which(motion_fd>max_fd);
pheno=pheno[-idx_fd,]
motion=motion[-idx_fd,]

motion_max=as.numeric(motion[,3])
motion_fd=as.numeric(motion[,33])
subject_ids=pheno[,1]

list_of_mats=list()
array_of_subjects=array()
count=0;



final_array_of_mc=array()
list_of_time_mats=list()
count=0;
    for (i in (1:length(subject_ids))) { ### iterate over subjects
        time_filename=paste('/scr/melisse1/NKI_enhanced/time_matrix_craddock200/0',subject_ids[i],'/mat.txt',sep='') # these are  correlation matrices from the coraddock 200 parcellation
        if (file.exists(file=time_filename)){
            z <- read.table(time_filename,sep=" ",dec=".",header=F)
            mat=do.call(cbind,z);
            if(length(which(is.nan(mat)==TRUE))<1){ ### no nan's
                if(length(which(mat==0))<1){ ### no 100 values. masked areas  with constant signal over time
                    idx=which(is.nan(mat)==TRUE)
                    mat[idx]=0;
                    count=count+1
                    idx=which(is.infinite(mat))
                    list_of_time_mats[[count]]=mat
                    array_of_subjects[count]=i;
                    final_array_of_mc=rbind(final_array_of_mc,motion[i,])
                    print(count)
                } else print(subject_ids[i])
            }
        }    else print(time_filename)
    }
final_array_of_mc=final_array_of_mc[-1,];
pheno=pheno[array_of_subjects,];
subject_ids=subject_ids[array_of_subjects];

source("/scr/melisse1/NKI/weight_conns.R")
window_width=77
bla_average=get_conn_mat(list_of_time_mats,mat)
bla=get_switch_mat(list_of_time_mats,mat,window_width)


#######cluster average connn
group_mean_average=apply(bla_average,c(2,3),mean)
group_mean_average1=apply(bla_average[seq(1,108,2),,],c(2,3),mean) #get mean, 108 is the number of subjects, this code is very unclean
group_mean_average2=apply(bla_average[seq(2,108,2),,],c(2,3),mean)
###thr95
group_mean_average_zeroed=group_mean_average
group_mean_average_zeroed[lower.tri(group_mean_average_zeroed,diag=T)]=0;
thr95=quantile(group_mean_average_zeroed,0.95,na.rm=T);
idx=which(group_mean_average_zeroed>thr95,arr.ind=T);
clust95=getLinkCommunities(cbind(idx,group_mean_average_zeroed[idx]),hcmethod='mcquitty')


#cex.main=3
#tiff("/scr/melisse1/NKI_enhanced/Colored_Dendo.tiff", width = 25, height= 15, units = 'in', res=150, compression = 'lzw')
#testa=plot(clust95,"summary")
#testa$hclust$order
#dev.off()
#for (i in (1:33)){
#a=1997
#if(length(which(testa$hclust$order[[a]]==clust95$clusters[[i]]))>0){
#    print(i)
#    print(a)
#    print(a+length(clust95$clusters[[i]]))
#}
#}
#tiff("/scr/melisse1/NKI_enhanced/Cluster_Related.tiff", width = 13, height= 9, units = 'in', res=300, compression = 'lzw')
#getClusterRelatedness(clust95,hcmethod='mcquitty',plot=TRUE)
#dev.off()

###hubanalyses
hist_hubs=hist(as.numeric(clust95$edgelist),breaks=201)
hubs <- rep(0, 200)
hubs[1:41]=hist_hubs$counts[1:41]
hubs[43:200]=hist_hubs$counts[42:199]
a=sort(as.numeric(names(clust95$numclusters)),index.return=T);
multihubs=clust95$numclusters[a$ix]
cor.test(hubs,multihubs)
write.csv(cbind(hubs,multihubs),'hubs.csv',row.names=F,col.names=F)

library("Rfit")
region_reg_sw_m=colMeans(region_switch_m)
lm1<-rfit(multihubs~hubs)
lm2<-rfit(region_reg_sw_m~hubs)

lm3<-rfit(hubs~multihubs)
lm4<-rfit(region_reg_sw_m~multihubs)
write.csv(cbind(hubs,multihubs,region_reg_sw_m,lm1$residuals,lm2$residuals,lm3$residuals,lm4$residuals),'/scr/melisse1/NKI_enhanced/hubs_dynamics.csv',row.names=F,col.names=F)

library("Rniftilib")# this is for plotting degree pictures`
parcell=nifti.image.read("/scr/melisse1/NKI_enhanced/correlation/craddock200_mni.nii.gz")
for (i in (1:200)){ # 200 ii the number of parcels
    idx=which(parcell[,,]==i,arr.ind=T)
    for (j in (1:(length(idx)/3))){
            parcell[idx[j,1],idx[j,2],idx[j,3]]=hubs[i]
        }
}
nifti.set.filenames(parcell, "/scr/melisse1/NKI_enhanced/normal_hubs.nii.gz")
nifti.image.write(parcell)

parcell=nifti.image.read("/scr/melisse1/NKI_enhanced/correlation/craddock200_mni.nii.gz")
for (i in (1:200)){
    idx=which(parcell[,,]==i,arr.ind=T)
        for (j in (1:(length(idx)/3))){
            parcell[idx[j,1],idx[j,2],idx[j,3]]=multihubs[i]
        }
}
nifti.set.filenames(parcell, "/scr/melisse1/NKI_enhanced/multi_hubs.nii.gz")
nifti.image.write(parcell)


###split half
idx_y1=which(as.numeric(pheno[,2])<25)[1:2:22]
idx_y2=which(as.numeric(pheno[,2])<25)[2:2:22]
idx_o1=which(as.numeric(pheno[,2])>45)[1:2:22]
idx_o2=which(as.numeric(pheno[,2])>45)[2:2:22]
group_mean_average1=apply(bla_average[idx_y1,,],c(2,3),mean)
group_mean_average2=apply(bla_average[idx_y2,,],c(2,3),mean)
group_mean_average3=apply(bla_average[idx_o1,,],c(2,3),mean)
group_mean_average4=apply(bla_average[idx_o2,,],c(2,3),mean)

clust95_1=getLinkCommunities(cbind(idx,group_mean_average1[idx]),hcmethod='mcquitty')
clust95_2=getLinkCommunities(cbind(idx,group_mean_average2[idx]),hcmethod='mcquitty')
clust95_3=getLinkCommunities(cbind(idx,group_mean_average3[idx]),hcmethod='mcquitty')
clust95_4=getLinkCommunities(cbind(idx,group_mean_average4[idx]),hcmethod='mcquitty')

coph_1=cophenetic(clust95_1$hclust);
coph_2=cophenetic(clust95_2$hclust);
coph_3=cophenetic(clust95_3$hclust);
coph_4=cophenetic(clust95_4$hclust);
mantel(coph_1,coph_2,permutations=20)
mantel(coph_1,coph_3,permutations=20)
mantel(coph_1,coph_4,permutations=20)

mantel(coph_4,coph_2,permutations=20)
mantel(coph_4,coph_2,permutations=20)
mantel(coph_4,coph_3,permutations=20)


group_mean_average1=apply(bla_average[seq(1,108,2),,],c(2,3),mean)
group_mean_average2=apply(bla_average[seq(2,108,2),,],c(2,3),mean)
clust95_1=getLinkCommunities(cbind(idx,group_mean_average1[idx]),hcmethod='mcquitty')
clust95_2=getLinkCommunities(cbind(idx,group_mean_average2[idx]),hcmethod='mcquitty')
coph_1=cophenetic(clust95_1$hclust);
coph_2=cophenetic(clust95_2$hclust);
library(vegan)
mantel(coph_1,coph_2,permutations=10)

A=true_edges(clust95);
#writeMat("/home/raid1/aschaefer/Matlab/Data/linkcomm/clust95.mat",A=A)

source("/scr/melisse1/NKI/weight_conns.R")
region_switch=array(0,dim= c(dim(bla)[1] ,dim(bla)[2] ))
region_switch=get_switching(clust95,bla)
region_switch_m=array(0,dim= c(dim(bla)[1] ,dim(bla)[2] ))
region_switch_m=get_multi_switching(clust95,bla)

subject_switching_m=rowMeans(region_switch_m);
subject_switching=rowMeans(region_switch);

cor.test(subject_switching_m,as.numeric(pheno[,2]),use='pairwise.complete.obs')

gender=array(0,dim(pheno)[1])
gender[which(pheno[,3]=="M")]=1

## plot regions included in analysis
pdf("ConnectivityPlot.pdf",width=7,height=7)
test=colMeans(region_switch)
plot(hubs,multihubs,title="Connectivity Plot",xlab="Connections",ylab="Networks", col = ifelse(test > 0.01,'red','black'), pch= ifelse(test > 0.01,20,18) )
legend(30,1, c("included in analysis","excluded from analysis"), col=c("red","black"),pch=c(20,18),cex=0.8);
dev.off()


#### plot dynamics to nifti
regional=colMeans(region_switch);
parcell=nifti.image.read("/scr/melisse1/NKI_enhanced/correlation/craddock200_mni.nii.gz")
for (i in (1:200)){
    idx=which(parcell[,,]==i,arr.ind=T)
        for (j in (1:(length(idx)/3))){
            parcell[idx[j,1],idx[j,2],idx[j,3]]=regional[i]*10
        }
}
nifti.set.filenames(parcell, "/scr/melisse1/NKI_enhanced/dynamic_hubs.nii.gz")
nifti.image.write(parcell)

regional_multi_switch=colMeans(region_switch_m);
parcell=nifti.image.read("/scr/melisse1/NKI_enhanced/correlation/craddock200_mni.nii.gz")
for (i in (1:200)){
    idx=which(parcell[,,]==i,arr.ind=T)
        for (j in (1:(length(idx)/3))){
            parcell[idx[j,1],idx[j,2],idx[j,3]]=regional_multi_switch[i]*10
        }
}
nifti.set.filenames(parcell, "/scr/melisse1/NKI_enhanced/dynamic_hubs_multi.nii.gz")
nifti.image.write(parcell)

###nodesize analysis
sizes = read.table("/scr/melisse1/Parcellations/craddock/nodesizes.txt", sep="\n")

source("/scr/melisse1/NKI_enhanced/pcor.R")


png(filename="/home/raid1/aschaefer/linkcomm/frontiers_motion.png", width=600,height=600)
plot()
plot(subject_switching_m,as.numeric(final_array_of_mc[,33]),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex.sub=1.5,cex=1.5,ylab="micro motion", xlab="switching")
res<-lm(as.numeric(final_array_of_mc[,33])~subject_switching_m)
abline(res)
write.csv(subject_switching_m,as.numeric(final_array_of_mc[,33]),"all_regions_motion.csv",row.names=F,col.names=F)

dev.off()

numcluster=clust95$numbers[3]
for (i in (seq(numcluster,numcluster*numcluster,numcluster))){
print(i/numcluster)
print(pcor.test(rowMeans(region_switch_m[,(i-(numcluster-1)):i]),rowMeans(mindw_new[,5:35]),cbind(as.numeric(pheno[,2]),as.numeric(final_array_of_mc[,33]))))
}

for (i in (seq(numcluster,numcluster*numcluster,numcluster))){
print(pcor.test(rowMeans(region_switch_m[,(i-(numcluster-1)):i]),vine,cbind(as.numeric(final_array_of_mc[,33]),as.numeric(pheno[,2]))))
}

motionmodel=lm(subject_switching~cbind(rowMeans(final_array_of_mc)));
cor.test(motionmodel$residuals,array_of_age,use='pairwise.complete.obs')


data = read.csv("/scr/melisse1/NKI_enhanced/pheno//mriq_all_subs_18_60.csv")


first_labels = c("Q01 - I thought about things I am currently worried about", 
    "Q02 - I thought about people I have just recently met", "Q03 - I thought of people I have known for a long time (friends)", 
    "Q04 - I thought about members of my family", "Q05 - I thought about an event that took place earlier today", 
    "Q06 - I thought about an interaction I may possibly have in the future", 
    "Q07 - I thought about an interaction with somebody that took place in the past", 
    "Q08 - I thought about something that happened at a place very close to me", 
    "Q09 - I thought about something that made me feel guilty", "Q11 - I thought about something that happened in the recent past (last couple of days but not today)", 
    "Q12 - I thought about something that happened a long time ago in the past", 
    "Q13 - I thought about something that made me angry", "Q14 - I thought about something that made me happy", 
    "Q15 - I thought about something that made me cheerful", "Q16 - I thought about something that made me calm", 
    "Q17 - I thought about something that made me sad", "Q18 - I thought about something that is important to me", 
    "Q19 - I thought about something that could still happen today", "Q20 - I thought about something that may take place in the distant future", 
    "Q21 - I thought about something that could take place in the near future (days or weeks but not today)", 
    "Q22 - I thought about personal worries", "Q23 - I thought about something that happened in a place far away from where I am now")

second_labels = c("Q24 - In the form of images:", "Q25 - In the form of words:", 
    "Q26 - Like an inner monologue or audiobook:", "Q27 - Like a television program or film:", 
    "Q28 - Had a strong and consistent personal narrative:", "Q29 - Had a clear sense of purpose:", 
    "Q30 - Vague and non-specific:", "Q31 - Fragmented and disjointed:)")
first_ind = c(4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 
    22, 23, 24, 25, 26)
second_ind = c(27, 28, 29, 30, 31, 32, 33, 34)
summary(data)

library(corrgram)
#corrgram(data, order = FALSE, lower.panel = panel.ellipse, upper.panel = panel.pie, text.panel = panel.txt, main = "Questionnaire relations")

fa.parallel(data[first_ind], fm = "pa", n.iter = 100)
first_fit = fa(data[,first_ind], nfactors = 5, scores = "tenBerge", fm = "pa", rotate = "oblimin")
fa.parallel(data[second_ind], fm = "pa", n.iter = 500)
second_fit = fa(data[,second_ind], nfactors = 3, scores = "tenBerge", fm = "pa", rotate = "oblimin")
regressors = data[c(1, 2, 3)]
regressors$firstSum = rowSums(data[first_ind])
regressors$secondSum = rowSums(data[second_ind])
regressors$allSum = regressors$firstSum + regressors$secondSum
first_scores = first_fit$scores
regressors = cbind(regressors, first_scores)
second_scores = second_fit$scores

regressors = cbind(regressors, second_scores)
#corrgram(regressors[c(-1)], order = FALSE, lower.panel = panel.ellipse, upper.panel = panel.pie,text.panel = panel.txt, main = "Questionnaire relations")
regressors$firstSum = rowSums(data[,first_ind])
regressors$secondSum = rowSums(data[,second_ind])
regressors$allSum = regressors$firstSum + regressors$secondSum
first_scores = first_fit$scores

regressor_new=matrix(nrow=dim(final_array_of_mc)[1],ncol=dim(regressors)[2])
count=1;
for (i in (1:length(subject_ids))){
    idx=which(as.numeric(data[,1])==as.numeric(subject_ids[i]));
    if(length(idx)>0){
        print(count)
        regressor_new[count,]=as.numeric(regressors[idx,]);
        count=count+1;
    }
}
idx=complete.cases(mindw_new[,5:35])
cor.test(subject_switching_m,regressor_new[,2])
for (i in (2:14)){
print(pcor.test(subject_switching_m,regressor_new[,i],cbind(as.numeric(final_array_of_mc[,33]),regressor_new[,2:3]), method = c("pearson")))
print(i)
}
source("/scr/melisse1/NKI_enhanced/pcor.R")
print(pcor.test(subject_switching_m,regressor_new[,8],cbind(as.numeric(final_array_of_mc[,33]),regressor_new[,2:3],regressor_new[,7],regressor_new[,9:14]), method = c("spearman")))
print(pcor.test(subject_switching_m,as.numeric(pheno[,2]),cbind(as.numeric(final_array_of_mc[,33])), method = c("spearman")))
print(pcor.test(region_switch[,84],regressor_new[,7],cbind(as.numeric(final_array_of_mc[,33]),regressor_new[,2:3],regressor_new[,8:14]), method = c("spearman")))

library("Rfit")
png(filename="/home/raid1/aschaefer/linkcomm/frontiers_positive.png", width=600,height=600)
lm1<-lm(subject_switching_m~as.numeric(final_array_of_mc[,33])+regressor_new[,2:3]+regressor_new[,7]+regressor_new[,9:14])
lm2<-lm(regressor_new[,8]~as.numeric(final_array_of_mc[,33])+regressor_new[,2:3]+regressor_new[,7]+regressor_new[,9:14])
plot(lm1$residuals,lm2$residuals,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex.sub=1.5,cex=1.5,ylab="positive mindwandering corrected", xlab="switching corrected")
res<-lm(lm2$residuals~lm1$residuals)
abline(res)
dev.off()
lm1<-rfit(subject_switching_m~as.numeric(final_array_of_mc[,33])+regressor_new[,2:3]+regressor_new[,7]+regressor_new[,9:14])
lm2<-rfit(regressor_new[,8]~as.numeric(final_array_of_mc[,33])+regressor_new[,2:3]+regressor_new[,7]+regressor_new[,9:14])
write.csv(cbind(lm1$residuals,lm2$residuals),"all_regions_positive.csv",row.names=F,col.names=F)

png(filename="/home/raid1/aschaefer/linkcomm/frontiers_age.png", width=600,height=600)
lm1<-rfit(subject_switching_m~as.numeric(final_array_of_mc[,33]))
lm2<-rfit(as.numeric(pheno[,2])~as.numeric(final_array_of_mc[,33]))
write.csv(cbind(lm1$residuals,lm2$residuals),"all_regions_age.csv",row.names=F,col.names=F)

write.csv(cbind(as.numeric(pheno[,2]),regressor_new[,8]),"age_sgt.csv",row.names=F,col.names=F)
lm1<-rfit(subject_switching_m~as.numeric(final_array_of_mc[,33])+regressor_new[,8])
lm2<-rfit(as.numeric(pheno[,2])~as.numeric(final_array_of_mc[,33])+regressor_new[,8])
write.csv(cbind(lm1$residuals,lm2$residuals),"all_regions_age_sgt.csv",row.names=F,col.names=F)



sw_l=array(0,108)
k=19
for (i in (1:108)){
sw_l[i]=sw_l[i]+sum(region_switch_m[i,k*32+(12)])+sum(region_switch_m[i,(12)*32+k])
}
for (i in (2:14)){
print(pcor.test(sw_l,regressor_new[,i],cbind(as.numeric(final_array_of_mc[,33]),regressor_new[,3],regressor_new[,2]), method = c("pearson")))
print(i)
}

second_category=regressor_new[,7:14]

local_p=array(200)
local_r=array(200)
for (i in (1:200)){
    local_p[i]=pcor.test(region_switch[,i],second_category[,7],cbind(as.numeric(final_array_of_mc[,33]),regressor_new[,2:3],second_category[,-7]), method = c("spearman"))$p.value
    local_r[i]=pcor.test(region_switch[,i],second_category[,7],cbind(as.numeric(final_array_of_mc[,33]),regressor_new[,2:3],second_category[,-7]), method = c("spearman"))$estimate
}
which(local_p<0.05/1600)

for (i in (1:14)){ #14 categories, age motion gender mindwandering and so on
    print(pcor.test(region_switch[,84],regressor_new[,i],cbind(as.numeric(final_array_of_mc[,33]),regressor_new[,2],regressor_new[,3]), method = c("pearson")))
    print(i)
}
print(pcor.test(region_switch[,84],regressor_new[,7],cbind(as.numeric(final_array_of_mc[,33]),regressor_new[,2:3],regressor_new[,8:14]), method = c("spearman")))

png(filename="/home/raid1/aschaefer/linkcomm/frontiers_past_84.png", width=600,height=600)
lm1<-lm(region_switch[,84]~as.numeric(final_array_of_mc[,33])+regressor_new[,3]+regressor_new[,2]+regressor_new[,8:14])
lm2<-lm(regressor_new[,7]~as.numeric(final_array_of_mc[,33])+regressor_new[,3]+regressor_new[,2]+regressor_new[,8:14])
plot(lm1$residuals,lm2$residuals,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex.sub=1.5,cex=1.5,ylab="mindwandering about the past corrected", xlab="switching corrected")
res<-lm(lm2$residuals~lm1$residuals)
abline(res)
dev.off()



lm1<-rfit(region_switch[,84]~as.numeric(final_array_of_mc[,33])+regressor_new[,2:3]+regressor_new[,8:14])# 8 to 14 are the 6 mindwandering categories
lm2<-rfit(regressor_new[,7]~as.numeric(final_array_of_mc[,33])+regressor_new[,2:3]+regressor_new[,8:14])
write.csv(cbind(lm1$residuals,lm2$residuals),"region84_past.csv",row.names=F,col.names=F)
cor(lm1$residuals,lm2$residuals,method='spearman')
print(pcor.test(region_switch[,84],regressor_new[,7],cbind(as.numeric(final_array_of_mc[,33]),regressor_new[,2:3],regressor_new[,8:14]), method = c("spearman")))

png(filename="/home/raid1/aschaefer/linkcomm/frontiers_visual_3.png", width=600,height=600)
lm1<-lm(region_switch[,3]~as.numeric(final_array_of_mc[,33])+regressor_new[,3]+regressor_new[,2])
lm2<-lm(regressor_new[,14]~as.numeric(final_array_of_mc[,33])+regressor_new[,3]+regressor_new[,2])
plot(lm1$residuals,lm2$residuals,cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex.sub=1.5,cex=1.5,ylab="mindwandering in images corrected", xlab="switching corrected")
res<-lm(lm2$residuals~lm1$residuals)
abline(res)
dev.off()

local_r=array(200)
for (i in (1:200)){
    local_r[i]=pcor.test(region_switch[,i],regressor_new[,14],cbind(as.numeric(final_array_of_mc[,33]),regressor_new[,2]), method = c("pearson"))$estimate
}
library("Rniftilib")
parcell=nifti.image.read("/scr/melisse1/NKI_enhanced/correlation/craddock200_mni.nii.gz")
for (i in (1:200)){
    idx=which(parcell[,,]==i,arr.ind=T)
    if( !(is.na(local_r[i]))) {
        if ( (length(idx>0)) ) {
            for (j in (1:(length(idx)/3))){
        parcell[idx[j,1],idx[j,2],idx[j,3]]=as.integer(local_r[i]*1000)
            }
        }
    }
}
nifti.set.filenames(parcell, "/scr/melisse1/NKI_enhanced/correlation/categ_14_corr_.nii.gz")
nifti.image.write(parcell)

cor.test(subject_switching_m[complete.cases(mindw_new[,5:35])],second_scores[,2])
pcor.test(subject_switching_m[idx],second_scores[,2],cbind(as.numeric(final_array_of_mc[idx,33])), method = c("kendall"))


