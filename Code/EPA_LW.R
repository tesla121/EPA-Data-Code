l2=list()
temp <<- list()
level=1
score<-function(Du1,Di,l,level,temp){
		if(level>(radius(sg))){
			return(0)
		}

	
	for(item in l){
		l2=append(l2,attributes(which(Ai[item,]==1))$names)

	}
	l2=unlist(l2)
	l3=setdiff(l2,l)
	l3=setdiff(l3,intersect(l3,temp))
	l3=as.character(l3)
	scr=0
	if(length(l3)==0){
		return (0)
	}
	for(item in l){
		scr=scr+(Di[item,item]/Du1[item,item])/(1/(1+log(Du1[item,item]))) 
	alpha=0
	temp=unique(unlist(append(temp,l)))
	if(length(l3)==0){
		return (0)
	}
	else{
		return (scr+score(Du1,Di,l3,level+1,temp))
	}
}



#Give paths to input underlying graphs and infection graphs
#Set the working directory
load(file="./EPA-Data-Code/Graphs/Facebook.RData")
load(file="./EPA-Data-Code/Infection Graphs (Single Source)/Facebook_Hetero_2.RData") 
#Facebook_Hetero_2 contains infected nodes list corresponding to the underlying graph, which in this case is Facebook. 
#Replace Facebook_Hetero_2 in the rest of the code according to the graph and infection size. For example if the 
#underlying graph is Regular and infection size is 40-60%, replace it with Regular_Hetero_40.

library(igraph)

#Code for Exoneration and Prominence based Age - LightWeight (EPA-LW) 



Au=as.matrix(get.adjacency(graph))
Du=diag(rowSums(Au))
rownames(Du)=as.character(V(graph)$name)
colnames(Du)=as.character(V(graph)$name)

sys=0
est_sum=0
EPA_LW_FB_Ht_2=list()
i=1
while(i<=100){
	source=V(graph)[Facebook_Hetero_2[[i]][1]]$name
	sg=induced_subgraph(graph, Facebook_Hetero_2[[i]], impl = c("copy_and_delete"))
	radius=radius(sg)
	nnodes<<-length(V(sg))
	Ai<<-as.matrix(get.adjacency(sg))
	Di=diag(rowSums(Ai))
	colnames(Di)=as.character(V(sg)$name)
	rownames(Di)=as.character(V(sg)$name)
	Du1=Du[rownames(Di),colnames(Di)]
	ecc=eccentricity(sg, vids = V(sg), mode = c("all"))
	ecc=degree(sg)
	estimated_ecc=V(sg)[as.numeric(which(min(ecc)==ecc))]$name
	estimated_ecc=V(sg)[as.numeric(which(max(ecc)==ecc))]$name
	scr_list=list()
	k=1
	for(node in as.character(estimated_ecc)) {
		scr=score(Du1,Di,as.character(node),level,temp)
		k=k+1
	} 
	scr_list=unlist(scr_list)
	index=which(max(scr_list)==scr_list)
	estimated=estimated_ecc[index]
	d=distances(sg, v = as.character(estimated), to = as.character(source), mode = c("all"), weights=NULL, algorithm = c("unweighted"))
	est_sum=est_sum+as.numeric(d[1])
	EPA_LW_FB_Ht_2[[i]]=as.numeric(d[1])
	if(i%%5==0){
		#save(EPA_LW_FB_Ht_2,file="./EPA-Data-Code/Result Objects/EPA_LW_FB_Ht_2.RData")
	
	}

	print(i)
	print(d)
	print(est_sum)
	i=i+1
	
}
