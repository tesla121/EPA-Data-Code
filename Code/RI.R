#Give paths to input underlying graphs and infection graphs
#Set the working directory
load(file="./EPA-Data-Code/Graphs/Facebook.RData")
load(file="./EPA-Data-Code/Infection Graphs (Single Source)/Facebook_Hetero_2.RData") 
#Facebook_Hetero_2 contains infected nodes list corresponding to the underlying graph, which in this case is Facebook. 
#Replace Facebook_Hetero_2 in the rest of the code according to the graph and infection size. For example if the 
#underlying graph is Regular and infection size is 40-60%, replace it with Regular_Hetero_40.

library(igraph)

#Code for Reverse Infection (RI)

RI_FB_Ht_2=list()
JC_sum=0
i=1
while(i<=100){
	source=V(graph)[Facebook_Hetero_2[[i]][1]]$name
	sg=induced_subgraph(graph, Facebook_Hetero_2[[i]], impl = c("copy_and_delete"))

id=list()
id2=list()

k=1
while(k<=length(V(sg))){

	id[[k]]=V(sg)[k]
	k=k+1
}

id2=list()
id2=id
stop=0

while(stop==0){
k=1
	while(k<=length(V(sg))) {
	neighbors=list()
	neighbors=unlist(adjacent_vertices(sg,k))
	for(n in neighbors){
		#print(n)
		id2[[n]]=unique(append(id2[[n]],id[[k]]))

	}
	k=k+1

}
id=id2


id=unique(id2)
#print(id)
j=1
while(j<=length(V(sg))) {
if(length(id[[j]])==length(V(sg))){
	stop=1
}
j=j+1
}

}

JCS=which(max(lengths(id))==lengths(id))
clo=closeness(sg, vids = JCS, mode = c( "all"),weights = NULL, normalized = TRUE)
index=which(max(clo)==clo)

	estimated_JC=V(sg)[as.numeric(JCS[index[1]])]$name

	d_JC=distances(sg, v = as.character(source), to = as.character(estimated_JC), mode = c("all"), weights=NULL, algorithm = c("unweighted"))
	print((d_JC))
	print(i)
	JC_sum=JC_sum+as.numeric(d_JC[1])
	RI_FB_Ht_2[[i]]=d_JC
	if(i%%5==0){
		#save(RI_FB_Ht_2, file="./EPA-Data-Code/Result Objects/RI_FB_Ht_2.RData")
	}
	i=i+1
}