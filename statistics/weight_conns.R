plot_hub_switch<-function(number,clust95,bla,subject) { #weights
    library(Hmisc)
    edges=true_edges(clust95)
    list_1=edges[(edges[,1]==number),]
    list_2=edges[(edges[,2]==number),]
    edges=rbind(list_1,list_2)
    print(edges)
    postscript(paste("/home/raid1/aschaefer/linkcomm/hub_switch_region_",number,"_subject",subject,".eps",sep=""), width=600,height=600)
    for (i in (unique(edges[,3]))) { ##loop through the clusters
        if(length(which(edges[,3]==i))>1){
            print(which(edges[,3]==i))
            print(paste("cluster:",i))
            print(length(which(edges[,3]==i)))
            z=array(0, dim=c(length(which(edges[,3]==i)),dim(bla)[4]))
            for (j in (1 : length(edges[which(edges[,3]==i),1]))){ ## loop through the connections
                print(j)
                #print(bla[1,edges[which(edges[,3]==i),1][j],edges[which(edges[,3]==i),2][j],])
                z[j,]=bla[subject,edges[which(edges[,3]==i),1][j],edges[which(edges[,3]==i),2][j],]
            }
            means=print(apply(z,2,mean))
            stdev=print(apply(z,2,sd)/j)
            par(ps=16)
            errbar(1:dim(bla)[4],means,means+stdev,means-stdev,col=i+9,errbar.col=i+9,ylim=c(-1.2,1.2),main="Switching Region",ylab="Average Connectivity (z)",xlab="Timepoint",cex=1.4,cex.axis=1.3,cex.lab=1.4)
            par(new=TRUE)
        }
    }
    dev.off()
}

write_hub_2_braingl<-function(number,clust95) { #weights
    edges=true_edges(clust95)
    list_1=edges[(edges[,1]==number),]
    list_2=edges[(edges[,2]==number),]
    edges=rbind(list_1,list_2)
    coord=read.table("/scr/melisse1/Parcellations/craddock/coordinates_200_mni_1mm.txt")
    print(unique(edges[,3]))
    print(edges)
    for (i in (unique(edges[,3]))) { ##loop through the clusters
        print(i)
        filename=paste("/scr/melisse1/NKI_enhanced/braingl/",number,"clust",i,".cxls",sep='')
        fileConn<-file(filename,"w")
        for (j in (1 : length(edges[which(edges[,3]==i),1]))){ ## loop through the connections
            print(j)
            line=(paste(c(coord[edges[which(edges[,3]==i),1][j],],coord[edges[which(edges[,3]==i),2][j],],1.0),collapse=' ') )
            writeLines(line,fileConn)
        }
        close(fileConn)           
    }
}

write_clust_2_braingl<-function(clust95) { #weights
    edges=true_edges(clust95)
    coord=read.table("/scr/melisse1/Parcellations/craddock/coordinates_200_mni_1mm.txt")
    for (i in (1 : max(edges[,3]))) { ##loop through the clusters
        filename=paste("/scr/melisse1/NKI_enhanced/braingl/clust",i,".cxls",sep='')
        fileConn<-file(filename,"w")
        for (j in (1 : length(edges[which(edges[,3]==i),1]))){ ## loop through the connections
            line=(paste(c(coord[edges[which(edges[,3]==i),1][j],],coord[edges[which(edges[,3]==i),2][j],],1.0),collapse=' ') )
            writeLines(line,fileConn)
        }
        close(fileConn)           
    }
}

weight_conns<-function(clust) { #weights 
    edges_m=data.matrix(clust$edgelist)
    clusters_m=data.matrix(clust$numclusters)
    weights=array(0,dim(length(edges_m)))
    for (i in (1:(dim(clust$edgelist)[1]))) {
            weights[i]=clusters_m[which(rownames(clusters_m)==edges_m[i,1])]+clusters_m[which(rownames(clusters_m)==edges_m[i,2])]
    }
    return(weights)
}


true_edges<-function(clust) { ##since there is a problem in the original r package with node relabeling when you have unconnected nodes, this version gives you a list of labeled edges with the "true" orginial node labels
    n_clusters=length(clust$clusters)
    true_edgelist=array(0,dim=c(dim(clust$edgelist)[1],3))
    counter=0;
    for (i in (1:(n_clusters))) {
            for (j in (1:(length(clust$clusters[[i]]) ))) {
                 edgenumber=clust$clusters[[i]][j];
                true_edgelist[edgenumber,1:2]=as.numeric(clust$edgelist[edgenumber,])
                true_edgelist[edgenumber,3]=i;
        }
    }
    return(true_edgelist)
}
get_conn_mat<-function(list_of_time_mats,mat) {
bla_average=array(0,dim=c(length(list_of_time_mats),dim(mat)[1],dim(mat)[1]))
for (j in (1:length(list_of_time_mats))){
    count=0;
    r=cor(t(list_of_time_mats[[j]]))
    average_z=0.5*(log(1+r) - log(1-r));
    idx=which(is.infinite(average_z)==TRUE);
    average_z[idx]=0;
	dim(average_z)
	dim(bla_average)
	bla_average[j,,]=average_z;
}
return(bla_average)
}

get_switch_mat<-function(list_of_time_mats,mat,window_width) {
bla=array(0,dim=c(length(list_of_time_mats),dim(mat)[1],dim(mat)[1],length(seq(window_width,dim(list_of_time_mats[[1]])[2],window_width))))
bla_average=array(0,dim=c(length(list_of_time_mats),dim(mat)[1],dim(mat)[1]))
for (j in (1:length(list_of_time_mats))){
    count=0;
    r=cor(t(list_of_time_mats[[j]]),use='complete')
    average_z=0.5*(log(1+r) - log(1-r));
    idx=which(is.infinite(average_z)==TRUE);
    average_z[idx]=0;
	print(dim(average_z))
	print(dim(bla_average))
	bla_average[j,,]=average_z;
    for (i in (seq(window_width,dim(list_of_time_mats[[1]])[2],window_width))) {### iteration over time
        count=count+1
        r=cor(t(list_of_time_mats[[j]][,(i-window_width+1):i]),use='complete');
        tmp_z=0.5*(log(1+r) - log(1-r));
        idx=which(is.infinite(tmp_z)==TRUE);
        tmp_z[idx]=0;
        bla[j,,,count]=average_z-tmp_z;
    }
}
return(bla)
}
get_switching<-function(clust,bla) { ##clust is an r linkcomm struct, bla is a 4 dimensional matrix with time changeing correlation matrices for all subjects, with dim(subjects, matrix dim1, matrix dim 2, time)
    cluster_list=true_edges(clust);
    region_switch=array(0,dim= c(dim(bla)[1] ,dim(bla)[2] ))
    for (s in (1:dim(bla)[1])) { ## run over subjects
        for (i in (1:dim(bla)[2])) { ## run through whole matrix number of regions
            idx1=which(cluster_list[,1]==i)
            idx2=which(cluster_list[,2]==i)
            list_of_ass_connections=rbind(cluster_list[idx2,],cluster_list[idx1,])
            occ=table(list_of_ass_connections[,3])
            if(length(occ)>1){##if there are two clusters
                secondmax=sort(occ,TRUE)[2];
                secondcluster=as.integer(names(secondmax))
                firstcluster=as.integer(names(sort(occ,TRUE)[1]))
                if((secondmax[[1]]>1)){#& (max(28==(c(firstcluster,secondcluster))))){ 
                     list_of_conn_idx1=which(list_of_ass_connections[,3]==firstcluster)
                    list_of_conn_idx2=which(list_of_ass_connections[,3]==secondcluster)
                    list_of_ass_connections[list_of_conn_idx1,1:2]
                    list_of_ass_connections[list_of_conn_idx2,1:2]
                    change1=array(0,dim=c(length(list_of_conn_idx1),dim(bla)[4]));
                    for (j in (1:length(list_of_conn_idx1))){  ##run trough associated edges (for the two largest clusters)
                        change1[j,]=bla[s,list_of_ass_connections[list_of_conn_idx1[j],1],list_of_ass_connections[list_of_conn_idx1[j],2],]
                    }
                    avg_change1=colMeans(change1)
                    change2=array(0,dim=c(length(list_of_conn_idx2),dim(bla)[4]));
                    for (j in (1:length(list_of_conn_idx2))){  ##run trough associated edges (for the two largest clusters)
                        change2[j,]=bla[s,list_of_ass_connections[list_of_conn_idx2[j],1],list_of_ass_connections[list_of_conn_idx2[j],2],]
                    }
                    avg_change2=colMeans(change2);
                    region_switch[s,i]=sum(abs(avg_change2-avg_change1));
                }
            }
        }
    }
return(region_switch)
}

get_multi_switching<-function(clust,bla) { ##clust is an r linkcomm struct, bla is a 4 dimensional matrix with time changeing correlation matrices for all subjects, with dim(subjects, matrix dim1, matrix dim 2, time)
    cluster_list=true_edges(clust);
    region_switch=array(0,dim= c(dim(bla)[1] ,clust$numbers[3]*dim(bla)[2] ))
    regional_switch=array(0,dim= c(dim(bla)[1] ,dim(bla)[2]))
    for (s in (1:dim(bla)[1])) { ## run over subjects
	print(s)
	count=0;
        for (i in (1:dim(bla)[2])) { ## run through whole matrix number of regions
            idx1=which(cluster_list[,1]==i)
            idx2=which(cluster_list[,2]==i)
            list_of_ass_connections=rbind(cluster_list[idx2,],cluster_list[idx1,])
            occ=table(list_of_ass_connections[,3])
            if(length(occ)>1){##if there are two clusters
            counter_cluster_comp=1;
            #print(sort(occ,TRUE))
        if((length(which(sort(occ,TRUE)>1))-1)>0){ ### not ideal but should work
		for (r in (1: (length(which(sort(occ,TRUE)>1))-1)) ){
            #print((length(which(sort(occ,TRUE)>1))))
            #print(r)
            #if((r+1)<length(which(sort(occ,TRUE)>1))){
			if((r+1)<=length(which(sort(occ,TRUE)>1))){ ### not ideal but should work
				for (k in ((r+1):length(which(sort(occ,TRUE)>1)))){
		        		secondmax=sort(occ,TRUE)[k];
					    secondcluster=as.integer(names(secondmax))
		        		firstcluster=as.integer(names(sort(occ,TRUE)[r]))
		        		if((secondmax[[1]]>1)){### if secondcluster is large enough
                            counter_cluster_comp=counter_cluster_comp+1;
					        list_of_conn_idx1=which(list_of_ass_connections[,3]==firstcluster)
					        list_of_conn_idx2=which(list_of_ass_connections[,3]==secondcluster)
					        list_of_ass_connections[list_of_conn_idx1,1:2]
					        list_of_ass_connections[list_of_conn_idx2,1:2]
					        change1=array(0,dim=c(length(list_of_conn_idx1),dim(bla)[4]));
					        for (j in (1:length(list_of_conn_idx1))){  ##run trough associated edges (for the two largest clusters)
						    change1[j,]=bla[s,list_of_ass_connections[list_of_conn_idx1[j],1],list_of_ass_connections[list_of_conn_idx1[j],2],]
					        }
					        avg_change1=colMeans(change1,na.rm=TRUE)
					        change2=array(0,dim=c(length(list_of_conn_idx2),dim(bla)[4]));
					        for (j in (1:length(list_of_conn_idx2))){  ##run trough associated edges (for the two largest clusters)
						    change2[j,]=bla[s,list_of_ass_connections[list_of_conn_idx2[j],1],list_of_ass_connections[list_of_conn_idx2[j],2],]
					        }
					        avg_change2=colMeans(change2,na.rm=TRUE);
						    count=count+1;
					        region_switch[s,firstcluster*(clust$numbers[3]-1)+secondcluster]=region_switch[s,firstcluster*(clust$numbers[3]-1)+secondcluster]+sum(abs(avg_change2-avg_change1));
                        regional_switch[s,i]=regional_switch[s,i]+sum(abs(avg_change2-avg_change1))*sqrt(length(list_of_conn_idx2)*length(list_of_conn_idx1))# account for number of connections, not over rate by sqrt
					    }
				}
            
			}
            
        }
        }
        regional_switch[s,i]=regional_switch[s,i]/(counter_cluster_comp) #account for number of cluster 
            
        }
        }
    }
return(regional_switch)
}
