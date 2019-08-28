#Give paths to input underlying graphs and infection graphs
#Set the working directory
load(file="./EPA-Data-Code/Graphs/Facebook.RData")
load(file="./EPA-Data-Code/Infection Graphs (Single Source)/Facebook_Hetero_2.RData") 
#Facebook_Hetero_2 contains infected nodes list corresponding to the underlying graph, which in this case is Facebook. 
#Replace Facebook_Hetero_2 in the rest of the code according to the graph and infection size. For example if the 
#underlying graph is Regular and infection size is 40-60%, replace it with Regular_Hetero_40.

library(igraph)

#Code for Dynamic Age (DA)

DA_FB_Ht_2=list() #Replace DA_FB_Ht_2 according to graph and infection size
DA_sum=0
i=1
while(i<=100){
	source=V(graph)[Facebook_Hetero_2[[i]][1]]$name
	sg=induced_subgraph(graph, Facebook_Hetero_2[[i]], impl = c("copy_and_delete"))
	A=as.matrix(get.adjacency(sg)) 

	l=list()
	j=1
	ev=max(eigen(A)$values)
	while(j<=nrow(A)){
		l[j]=abs(ev-max(eigen(A[-j,-j])$values))/ev
		j=j+1
		#print(j)
	}
	l=unlist(l)
	estimated_DA=V(sg)[as.numeric(which(max(l)==l))]$name

	d_DA=distances(sg, v = as.character(source), to = as.character(estimated_DA), mode = c("all"), weights=NULL, algorithm = c("unweighted"))
	print((d_DA))
	print(i)
	DA_sum=DA_sum+as.numeric(d_DA[1])
	DA_FB_Ht_2[[i]]=d_DA
		if(i%%5==0){
			#save(DA_FB_Ht_2, file="./EPA-Data-Code/Result Objects/DA_FB_Ht_2.RData")
		}
	i=i+1
}
