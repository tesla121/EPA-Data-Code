#Give paths to input underlying graphs and infection graphs
#Set the working directory
load(file="./EPA-Data-Code/Graphs/Facebook.RData")
load(file="./EPA-Data-Code/Infection Graphs (Single Source)/Facebook_Hetero_2.RData") 
#Facebook_Hetero_2 contains infected nodes list corresponding to the underlying graph, which in this case is Facebook. 
#Replace Facebook_Hetero_2 in the rest of the code according to the graph and infection size. For example if the 
#underlying graph is Regular and infection size is 40-60%, replace it with Regular_Hetero_40.

library(igraph)
library(RSpectra)

#Code for Minimum Description Length (MDL) - Fast approach

	A=as.matrix(get.adjacency(graph)) 
	D=diag(rowSums(A))
	L=D-A
	mdl_FB_Ht_2=list() #Replace mdl_FB_Ht_2 according to graph and infection size
	mdl_sum=0
	i=1
	while(i<=100){
		source=V(graph)[Facebook_Hetero_2[[i]][1]]$name
		sg=induced_subgraph(graph, Facebook_Hetero_2[[i]], impl = c("copy_and_delete"))
		Lm=L[c(Facebook_Hetero_2[[i]]),c(Facebook_Hetero_2[[i]])]
		e=eigs(Lm,1,which="SM")
		si=which(max(abs(e$vector))==abs(e$vector))		
		x=Facebook_Hetero_2[[i]][si]
		estimated_mdl=V(graph)[x]$name
		d_mdl=distances(sg, v = as.character(source), to = as.character(estimated_mdl), mode = c("all"), weights=NULL, algorithm = c("unweighted"))
		print(d_mdl)
		print(i)
		mdl_sum=mdl_sum+as.numeric(d_mdl[1])
		mdl_FB_Ht_2[[i]]=d_mdl
		if(i%%5==0){
			#save(mdl_FB_Ht_2, file="./EPA-Data-Code/Result Objects/mdl_FB_Ht_2.RData")
		}
		i=i+1
	}

