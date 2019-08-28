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
	#print(l3)
	l3=as.character(l3)
	scr=0
	if(length(l3)==0){
		return (0)
	}
	for(item in l){
		scr=scr+(Di[item,item]/Du1[item,item])/(1/(1+log(Du1[item,item]))) 
	}
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
load(file="./EPA-Data-Code/Graphs/Jazz.RData")
load(file="./EPA-Data-Code/Infection Graphs (Multiple Source)/2/Jazz_Hetero_10_2S.RData") 
#Jazz_Hetero_10_2S contains infected nodes list corresponding to the underlying graph, which in this case is Jazz. 

library(igraph)
library(gtools)

#Code for EPA_SSI


Au=as.matrix(get.adjacency(graph))
Du=diag(rowSums(Au))
rownames(Du)=as.character(V(graph)$name)
colnames(Du)=as.character(V(graph)$name)

nos=2 #Initialize the number of sources
EPA_SSI_Jazz_Ht_10_2S_DR=list()
EPA_SSI_Jazz_Ht_10_2S_AED=list()
f1_sum=0
out_min_sum=0
sys=0
est_sum=0
EPA_SSI_Jazz_Ht_10_2S_DR=list()
i=1
while(i<=100){
	sources=V(graph)[Jazz_Hetero_10_2S[[i]][1:nos]]$name
	sg=induced_subgraph(graph, Jazz_Hetero_10_2S[[i]], impl = c("copy_and_delete"))
	radius=radius(sg)
	nnodes<<-length(V(sg))
	Ai<<-as.matrix(get.adjacency(sg))
	Di=diag(rowSums(Ai))
	colnames(Di)=as.character(V(sg)$name)
	rownames(Di)=as.character(V(sg)$name)
	Du1=Du[rownames(Di),colnames(Di)]
	kk=1
	seeds=list()
	est_sum_in=0
	old=0
	while(kk<=nos){
		scr_list=list()
		k=1
		for(node in as.character(V(sg)$name) ){
			scr=score(Du1,Di,as.character(node),level,temp)
			penalty=as.numeric(eccentricity(sg, vids = node, mode = c("all")))
			scr_list[k]=scr/penalty
			k=k+1
		} 
		scr_list=unlist(scr_list)
		index=which(max(scr_list)==scr_list)[1]
		seeds=append(seeds,as.character(V(sg)[index]$name))
		seeds=unlist(seeds)
		Ai[which(rownames(Ai)==(V(sg)[index]$name)),]=c(rep(0,length(rownames(Ai))))
		Ai[,which(rownames(Ai)==(V(sg)[index]$name))]=c(rep(0,length(rownames(Ai))))
		Di[which(rownames(Di)==(V(sg)[index]$name)),]=c(rep(0,length(rownames(Di))))
		Di[,which(rownames(Di)==(V(sg)[index]$name))]=c(rep(0,length(rownames(Di))))
		Du1=Du[rownames(Di),colnames(Di)]
		if(length(Di)==0){
		#break
		}
		kk=kk+1
		new=max(scr_list)
		if(new==old){
		#break
		}
		else{
			old=new
		}
	}
	outer_sum=10000000000
perm=permutations(n=nos,r=nos,v=seeds)
pit=1   
while(pit<=nrow(perm))
{
	#print(perm[pit,])
	sit=1
	inner_sum=0
	while(sit<=length(perm[pit,])){
		inner_sum=inner_sum+as.numeric(distances(sg, v = as.character(sources[sit]), to = as.character(perm[pit,sit]), mode = c("all"), weights=NULL, algorithm = c("unweighted")))
		sit=sit+1
	}
	inner_sum=inner_sum/length(sources)
	#print(inner_sum)
	if(inner_sum<outer_sum){
		outer_sum=inner_sum
	}
	pit=pit+1

}
print(i)
print(outer_sum)
est_sum=est_sum+outer_sum

print(est_sum)

dr=length(intersect(as.character(sources),as.character(seeds)))/length(sources)
print(dr)
print("_______________")
EPA_SSI_Jazz_Ht_10_2S_DR[[i]]=dr
save(EPA_SSI_Jazz_Ht_10_2S_DR, file="./EPA-Data-Code/Result Objects/EPA_SSI_Jazz_Ht_10_2S_DR.RData")

EPA_SSI_Jazz_Ht_10_2S_AED[[i]]=outer_sum
save(EPA_SSI_Jazz_Ht_10_2S_AED, file="./EPA-Data-Code/Result Objects/EPA_SSI_Jazz_Ht_10_2S_AED.RData")
	i=i+1
	
}
