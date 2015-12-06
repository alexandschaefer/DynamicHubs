library(WGCNA)
library(psych)
library(R.matlab)
library(e1071)
library(calibrate)
library(gplot)
library(sciplot)
### this is an example wrapper file
###
### before you can use this code you have to 
### create textfiles containing timeseries data 
### the timeseries should be extracted using the 
### craddock 200 parcellation
###
### the code will create text files containing 
### regional switching information and averaged 
### switching information for every subject
###
### the switching variables cover the
### switching infomration between the 
### two largest networks
###
### in the publication we used the _m variables
### which contain the switching information
### between all networks



### read in the data
subject_ids<- as.matrix(read.csv(file="../text_file_with_filenames.txt")) #provide a textfile with the absolute filenames
list_of_time_mats={}
count=0;
    for (i in (1:length(subject_ids))) { ### iterate over subjects
        time_filename=paste('',subject_ids[i],sep='') # these are  correlation matrices from the coraddock 200 parcellation
        if (file.exists(file=time_filename)){
            z <- read.table(time_filename,sep=" ",dec=".",header=F)
            mat=do.call(cbind,z);
            if(length(which(is.nan(mat)==TRUE))<1){ ### no nan's
                    idx=which(is.nan(mat)==TRUE)
                    mat[idx]=0;
                    count=count+1
                    idx=which(is.infinite(mat))
                    list_of_time_mats[[count]]=mat
                    #array_of_subjects[count]=i;
                    print(count)
            }
        }    else print(time_filename)
    }
subject_ids=subject_ids[array_of_subjects];

###hubanalyses
load("statistics/publication_result/clust95.RData")
source("statistics/weight_conns.R")
window_width=20 ## set the timewindow, this depends on your TR. e.g. TR=2.3, window width 20 -> 46sec
bla_average=get_conn_mat(list_of_time_mats,mat) ## compute the average connectivity mat
bla=get_switch_mat(list_of_time_mats,mat,window_width) ### compute the connectiovity mats for all timepoints


region_switch=array(0,dim= c(dim(bla)[1] ,dim(bla)[2] ))
region_switch=get_switching(clust95,bla) ### compute regional switching
region_switch_m=array(0,dim= c(dim(bla)[1] ,dim(bla)[2] ))
region_switch_m=get_multi_switching(clust95,bla) ### compute multi switching

subject_switching_m=rowMeans(region_switch_m,na.rm=TRUE); ## this is just an average over the subject
subject_switching=rowMeans(region_switch,na.rm=TRUE); ### same as above


### write out the results 
write.csv(subject_switching,'subject_switching_allrois.csv',row.names=F,col.names=F)
write.csv(subject_switching_m,'subject_switching_m_allrois.csv.csv',row.names=F,col.names=F)
write.csv(region_switch,'region_switching.csv_allrois.csv',row.names=F,col.names=F)
write.csv(region_switch_m,'region_switching_m.csv_allrois.csv',row.names=F,col.names=F)




