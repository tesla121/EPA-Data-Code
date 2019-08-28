library(igraph)


load(file="./EPA-Data-Code/Graphs/Jazz.RData")



Jazz_Hetero_10_2S=list()

k=1
while(k<=100){ #Generate 100 infection graphs
	random=list()
	print(k)
	sample=sample(V(graph))
random=as.numeric(c(sample[1],sample[2])) #Randomly pick two sources

t=0
neighbors=list()
temp_active=list()
perm_active=list()
temp=list()
adj=list()
perm_active=random

while(1)
{
	t=t+1
	neighbors=list()
	for(node in unlist(perm_active))
	{
		adj=adjacent_vertices(graph,(node))
		i=1
		while(i<=length(unlist(adj)))
		{
			neighbors[length(unlist(neighbors))+1]=unlist(adj)[i]
			i=i+1
		}
	}
	neighbors=unique(unlist(neighbors))
	i=1	
	rel_neighbors=setdiff(unlist(neighbors),unlist(perm_active))
	#print(t)
	#print(neighbors)
	#print(rel_neighbors)
	temp=perm_active
	#temp_active=list()
	for(inactive in rel_neighbors)
	{	
		count=0
		edge_list=list()
		for(active in temp)
		{	
			
			
			if(are.connected(graph,(inactive),(active)))
			{
				count=count+1
			#	print("c")
				#print(count)
				e=get.edge.ids(graph, c((inactive),(active)))
				edge_list[length(unlist(edge_list))+1]=E(graph)[e]$weight

				#print(unlist(edge_list))
				#print(active)
				#print(inactive)
			}
			#print(edge_list)
		}
		if(count>1)
		{
			i=1
			prob=1
			while(i<=count)
			{
				prob=prob*(1-unlist(edge_list)[i])
			#	print(prob)
				i=i+1
			}
			prob=1-prob
			#print(prob)
		}
		else
		{

			prob=unlist(edge_list)
		}
		if(prob>runif(n=1, min=0, max=1))
		{
			perm_active[length(unlist(perm_active))+1]=(inactive)
			
			
			#print(unlist(perm_active))				
			#print(prob)
		}
	}
	#print(length(unlist(perm_active)))
	if(length(unlist(perm_active))>=length(V(graph))*0.10) #Set infection size
			{
				#if(length(unlist(perm_active))<=length(V(graph))*0.05){
						sg=induced_subgraph(graph, perm_active, impl = c("copy_and_delete"))
						if(is.connected(sg)==F) {
							print("NOT CONNECTED")
							break
						}
						else{
				Jazz_Hetero_10_2S[[length(Jazz_Hetero_10_2S)+1]]=unlist(perm_active)
				print(unlist(length(perm_active)))	
			#}	
			#save(Jazz_Hetero_10_2S, file="./EPA-Data-Code/Infection Graphs (Multiple Source)/2/Jazz_Hetero_10_2S.RData")

				k=k+1
				break
			}
			}		
}
}
