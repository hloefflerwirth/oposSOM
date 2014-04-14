


if( length(set.list$spots) >= 3 )
{
# 	
# 	
# 	###    signifikanzfilter
# 	
# 		significance_filter = T
# 	
# 	
# 		r=c()
# 		t=c()
# 		p=c()
# 		fdr.threshold=0.05
# 	
# 		for(i in 1: nrow(indata))
# 		{
# 			k=som.nodes[i]
# 			suppressWarnings({ r[i]=cor(indata[i,],metadata[k,]) })
# 			t[i]=r[i]/(sqrt((1-r[i]^2)/(ncol(indata)-2)))
# 			p[i]=1-pt(t[i], ncol(indata)-2)
# 		}
# 		names(p)=rownames(indata)
# 		insig=names(p)[which(is.na(p))]
# 		p1=p[which(!is.na(p))]
# 		suppressWarnings({  f=fdrtool(p1,statistic="pvalue",plot=F,verbose=F)  })
# 		insig=c(insig,names(p)[which(f$lfdr>fdr.threshold)])
# 	
# 		
# 	###           b_ij+a_ij
# 	###    w=-------------------
# 	###      min(k_i,k_j)+1-a_ij
# 	
# 	###    adj_matrix = a
# 	
# 	
# 		adj_matrix=cor(t(metadata),t(metadata))
# 		
# 		K=c()
# 	
# 		for (i in 1: nrow(metadata))
# 		{
# 			if(significance_filter)
# 			{
# 				ki=names(som.nodes[-match(insig,rownames(indata))])[which(som.nodes[-match(insig,rownames(indata))]==i)]
# 			}else
# 			{
# 				ki=names(som.nodes)[which(som.nodes==i)]
# 			}
# 			K[i]=length(ki)
# 		}
# 	
# 	
# 		if(significance_filter)
# 		{
# 			f=t(t(K))%*%t(K)
# 			f=f/max(f)
# 		}else
# 		{
# 			f=1
# 		}
# 	
# 
# 		diag(adj_matrix) = 0
# 		adj_matrix=adj_matrix*f
# 	
# 	
# 	
# 		w=matrix(NA,nrow=nrow(metadata),ncol=nrow(metadata))
# 	
# 		b = adj_matrix %*% adj_matrix
# 	
# 		k_sum = apply(adj_matrix,2,function(x){sum(abs(x))})
# 	
# 		for(i in 1: nrow(metadata))
# 		{
# 			for(j in 1: nrow(metadata))
# 			{
# 				w[i,j] = ( b[i,j] + adj_matrix[i,j] )   /   ( min(k_sum[i],k_sum[j]) + 1 - abs(adj_matrix[i,j]) )
# 			}
# 		}
# 	
# 	
# 	
# 	
# 	###    metagene groups
# 	
# 	###    M: groesse der spots
# 	
# 		M=c()
# 		for (i in 1:length(set.list$spots) )
# 		{
# 			M[i]=length(which(set.list$spots[[i]]$mask==1))
# 		}
# 	
# 		omega=matrix(NA,length(M),length(M))
# 		
# 		for( i in 1: length(M))
# 		{
# 			for( j in 1: length(M))
# 			{
# 				omega[i,j]=1/(M[i]*M[j])*sum(w[which(set.list$spots[[i]]$mask==1),which(set.list$spots[[j]]$mask==1)])
# 			}
# 		}
# 	
# 		rownames(omega)=names(set.list$spots)
# 		colnames(omega)=names(set.list$spots)
# 	
# 		
# 	###    normierung von omega auf Intervall -1..1
# 	
# 		if( max(omega[which(omega>0)])-min(omega[which(omega>0)]) )	
# 			omega[which(omega>0)]=(omega[which(omega>0)]-min(omega[which(omega>0)]))/(max(omega[which(omega>0)])-min(omega[which(omega>0)]))
# 		else
# 			omega[which(omega>0)]=1
# 		
# 		if( max(omega[which(omega<0)])-min(omega[which(omega<0)]) )	
# 			omega[which(omega<0)]=(omega[which(omega<0)]-min(omega[which(omega<0)]))/(max(omega[which(omega<0)])-min(omega[which(omega<0)]))-1
# 		else
# 			omega[which(omega<0)]=1
# 		
# 		
# 	###    threshold durch 25- und 75-percentile
# 	
# 		omega[intersect(which(omega> quantile(omega,c(0.15,0.85))[1]),which(omega< quantile(omega,c(0.15,0.85))[2]))]=0
# 	
# 	
# 	 	for(i in 1: ncol(omega))
# 	 	{
# 	 		for(j in 1: ncol(omega))
# 	 		{
# 	 			if(j>=i)
# 	 			{
# 	 				omega[i,j]=0
# 	 			}
# 	 		}
# 	 	}
# 	
# 		omega2=omega
# 		omega2[which(omega2!=0)]=1
# 	




		
		
		
		adj_matrix=cor(t(set.list$spotdata),t(set.list$spotdata))
		diag(adj_matrix) = 0
		
		
		omega = matrix(NA,nrow=nrow(set.list$spotdata),ncol=nrow(set.list$spotdata), dimnames=list(names(set.list$spots),names(set.list$spots)))	
		b = adj_matrix %*% adj_matrix
		k_sum = apply(adj_matrix,2,function(x){sum(abs(x))})
		
		for(i in 1: nrow(set.list$spotdata))
		{
			for(j in 1: nrow(set.list$spotdata))
			{
				omega[i,j] = ( b[i,j] + adj_matrix[i,j] )   /   ( min(k_sum[i],k_sum[j]) + 1 - abs(adj_matrix[i,j]) )
			}
		}
		diag(omega) = 0
		
		omega[ which( abs(omega) < 0.5 ) ] = 0
		





		g <- graph.adjacency( omega, weighted=T,  mode="undirected" )

		layout = layout.fruchterman.reingold(g)
		layout = layout.norm( layout, xmin=-1, xmax=1, ymin=-1, ymax=1 )

		V(g)$label = names(set.list$spots)

		n.spot = apply( sample.spot.matrix, 1, sum ) 
		V(g)$size = 8*(n.spot-min(n.spot))/(max(n.spot)-min(n.spot))+6


		
		if(	length(E(g)) > 0	)
		{	
			E(g)$color = c("red2","green2")[ 1.5 + sign( E(g)$weight ) / 2 ]
	
			E(g)$width = abs( E(g)$weight )
			E(g)$width = 5*(E(g)$width-min(E(g)$width))/(max(E(g)$width)-min(E(g)$width))+0.5
		}	
	
		
	
		par(mar=c(1,1,1,1), mfrow=c(1,1) )
		
		plot( g, layout=layout, vertex.label.color="black" ,vertex.label.cex=1, vertex.color="grey", main="", xlim=c(-1,1), ylim=c(-1.1,1.1) ) 
			legend("bottomright", c("positive correlation", "negative correlation"), lty = c(1,1), lwd = 4,col=c("green2","red2")) 
			text( -1.65, 1.1, "WTO network", adj=0, cex=2 )
		

		par(mar=c(4.9,13.5,4.9,13.5))
		image( matrix(set.list$overview.mask, preferences$dim.som1), col="gray90", axes=F )
		par(new=T, mar=c(1,1,1,1))
		plot( g, layout=t( sapply( set.list$spots, function(x) x$position ) ), vertex.label.color="black" ,vertex.label.cex=1, vertex.color="grey", main="", xlim=c(-1.2,1.2), ylim=c(-1.2,1.2) )		
			rect(-1.1,-1.1,1.1,1.1)
			legend("bottomright", c("positive correlation", "negative correlation"), lty = c(1,1), lwd = 4,col=c("green2","red2")) 

		
		


		if(	length(E(g)) > 0	)
		{	

			g2 = g - E(g)[ weight < 0 ]
		
			par(mar=c(4.9,13.5,4.9,13.5))
			image( matrix(set.list$overview.mask, preferences$dim.som1), col="gray90", axes=F )
			par(new=T, mar=c(1,1,1,1))
			plot( g2, layout=t( sapply( set.list$spots, function(x) x$position ) ), vertex.label.color="black" ,vertex.label.cex=1, vertex.color="grey", main="", xlim=c(-1.2,1.2), ylim=c(-1.2,1.2) )		
				rect(-1.1,-1.1,1.1,1.1)
				legend("bottomright", "positive correlation", lty = 1, lwd = 4,col="green2") 
			
	
	
			g2 = g - E(g)[ weight > 0 ]
	
			par(mar=c(4.9,13.5,4.9,13.5))
			image( matrix(set.list$overview.mask, preferences$dim.som1), col="gray90", axes=F )
			par(new=T, mar=c(1,1,1,1))
			plot( g2, layout=t( sapply( set.list$spots, function(x) x$position ) ), vertex.label.color="black" ,vertex.label.cex=1, vertex.color="grey", main="", xlim=c(-1.2,1.2), ylim=c(-1.2,1.2) )		
				rect(-1.1,-1.1,1.1,1.1)
				legend("bottomright", "negative correlation", lty = 1, lwd = 4,col="red2") 	
			
		}
		
		
}		




