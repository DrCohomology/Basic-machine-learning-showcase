library(readr)
library(imputeTS)
library(haven)
library(gtools)
library(Hmisc)
library(magrittr)
library(ggplot2)
library(ggcorrplot)
library(PerformanceAnalytics)
library(factoextra)
library(cluster)

# ---- Data Loading ----


data <- read.csv("C:/Users/Io/Desktop/SML Unsupervised/water_treatment.csv", sep=',', header=FALSE, stringsAsFactors=FALSE)
length(data[data=='?'])
data[data=='?'] <- NA
data <- data %>% subset(select=-c(1)) %>% data.frame

for(i in 1:ncol(data)){
	data[[i]] <- as.numeric(data[[i]])
	data_mean <- mean(data[[i]], na.rm=TRUE)
	data[[i]] <- data[[i]] %>% na.replace(data_mean) 
}
data <- scale(data) %>% data.frame()

boxplot(data, macro.labels='c')

p <- ncol(data)
n <- nrow(data)

# ---- Correlation, number of highly correlated pairs ----

cor_matrix <- cor(data)
length(cor_matrix[cor_matrix != 1 & (cor_matrix >= 0.8 | cor_matrix <= -0.8) ])/2

# chart.Correlation(data, pch=19, subset=cor_matrix[cor_matrix != 1 & (cor_matrix >= 0.8 | cor_matrix <= -0.8) ])

# ---- PCA ----

eig_val <- eigen(cor_matrix)$values
eig_vec <- eigen(cor_matrix)$vectors

eig_val_CI <- 0.95 * sqrt(2/n) * eig_val + eig_val

expl_var <- eig_val/p
CI_var <- eig_val_CI/p

cumexpl_var <- cumsum(expl_var)
cumCI_var <- cumsum(CI_var)
tab <- cbind(eig_val, expl_var, CI_var, cumexpl_var, cumCI_var) %>% round(3)
colnames(tab) <- c('eigenvalue', 'variance', 'variance_CI', 'cumulative_variance', 'cumulative_variance_CI')
tab

COL <- (eig_val_CI >= 1) 
COL[COL==TRUE] <- 'red'
COL[COL==FALSE] <- 'black'

plot(eig_val_CI, type='b', pch=19, main='scree plot', xlab='component', ylab='eigenvalue', col=COL)
abline(h=1, col="red")

# 10 components selected with updated Kaiser Rule, corresponding to 0.851 of total variance
# 13 components selected with fraction of explained variance > 0.9

P <- eig_vec[,1:10] 
PCA_data <- as.matrix(data) %*% P %>% scale

# ---- H Clustering ----

ds <- dist(PCA_data)

clusters <- function(H){
	data.frame( row.names = paste0('Cluster', seq_along(H$height)),
				height = H$height,
			   	components = ifelse(H$merge<0, H$labels[abs(H$merge)], paste0('Cluster', H$merge)),
				stringsAsFactors = FALSE
		)
}

sin_linkage <- hclust(ds, method='single')
#plot(sin_linkage, main='single')
#rect.hclust(sin_linkage, k=5)

cpl_linkage <- hclust(ds, method='complete')
#plot(cpl_linkage, main='complete')

# ---- Optimal Number of clusters ----

optimal_number <- function(ds, method, max_K, target='wss', plot=FALSE){
	
	valid_target <- list('within_variance', 'wss')
	
	out <- list()
	number_clust <- c()
	obj <- c()
	
	
	clust <- hclust(ds, method=method)
	numbers <- 2:max_K
	
	for(K in numbers){
		data$cluster <- cutree(clust, k=K)
		temp <- rep(NA, K)
			
		# minimize target 
		for(k in 1:K) {
			clusterk <- data[data$cluster == k, ]
			clusterk <- clusterk[, -c(length(clusterk))]	# drop the clusters column
			
			if(target=='within_variance'){
				if(nrow(clusterk) > 1)
					temp[k] <- sum((clusterk - colMeans(clusterk))^2) /(nrow(clusterk)-1)
				else
					temp[k] <- 0
			}
			if(target=='wss')
				temp[k] <- sum((clusterk - colMeans(clusterk))^2)
		}

		obj <- obj %>% append(sum(temp))
	}
	
	opt <- min(numbers[obj == min(as.numeric(obj))])
	
	if(plot){
		plot(numbers, type='l', obj, 
			xlab='number of clusters', ylab=target,
			main=paste0(method, ' linkage', '\n', 
						'optimal number of clusters: ', opt, '\n', 
						target, ': ', obj[numbers == opt] %>% round(2)
			)
		)
		abline(v=opt, col='red')
	}
	
	out$k <- opt
	
	if(target=='within_variance')
		out$summary <- data.frame('K' = numbers, 'within_variance' = obj)
	if(target=='wss')
		out$summary <- data.frame('K' = numbers, 'wss' = obj)
	
	return(out)
}

sin_opt <- optimal_number(ds, 'single', 200, 'wss', plot=FALSE)
cpl_opt <- optimal_number(ds, 'complete', 200, 'wss', plot=FALSE)

plot(sin_opt$summary$K, sin_opt$summary$wss, type='l', 
	main='single linkage', xlab='number of clusters', ylab='wss'
	); abline(v=30, col='red')
plot(cpl_opt$summary$K, cpl_opt$summary$wss, type='l', 
	main='complete linkage', xlab='number of clusters', ylab='wss'
	); abline(v=12, col='red')

sin_opt$k <- 30
cpl_opt$k <- 12

plot(sin_linkage, main=paste0('single linkage \n optimal number of clusters : ', sin_opt$k), labels=FALSE); 
rect.hclust(sin_linkage, k=sin_opt$k)
plot(cpl_linkage, main=paste0('complete linkage \n optimal number of clusters : ', cpl_opt$k), labels=FALSE); 
rect.hclust(cpl_linkage, k=cpl_opt$k)

# ---- Optimal dataset ----

sin_top <- cutree(sin_linkage, k=sin_opt$k)
cpl_top <- cutree(cpl_linkage, k=cpl_opt$k)

data_sin <- data.frame(PCA_data, 'cluster'=sin_top)
data_cpl <- data.frame(PCA_data, 'cluster'=cpl_top)

aggregate(data_sin$cluster, by=list(data_sin$cluster), FUN=length)
aggregate(data_cpl$cluster, by=list(data_cpl$cluster), FUN=length)

# ---- Cluster Visualization ----

data_single <- list('data'=data, 'cluster'=sin_top)
data_complete <- list('data'=data, 'cluster'=cpl_top)

fviz_cluster(data_single, 
			 main = 'single linkage',
			 frame.type = 'convex', 
			 geom = 'point', 
			 labelsize = 1,
) + theme_minimal() 

fviz_cluster(data_complete, 
			 main = 'complete linkage',
			 frame.type = 'convex', 
			 geom = 'point', 
			 labelsize = 1,
) + theme_minimal()

# ---- elbow kmeans----

n <- nrow(PCA_data)
times <- 20
TIMES <- 40
N <- n%/%10
wss <- rep(NA, N)
wss_std <- rep(NA, N)
wssK <- rep(NA, times)
numbers <- 1:N
elbow <- rep(NA, TIMES)


for(T in 1:TIMES){
	for(K in numbers){
		for(t in 1:times)
			wssK[t] <- kmeans(PCA_data, centers=K)$tot.withinss
		wss[K] <- mean(wssK)
		wss_std[K] <- sd(wssK)
	}
	
	for(x in 1:(length(wss)-1)){
		if(wss[x]-wss_std[x] < wss[x+1] + wss_std[x+1]){
			elbow[T] <- x
			break
		}
	}
}
boxplot(elbow, main=paste0('elbow \n mean: ', mean(elbow)))

qplot(numbers, wss, main=paste0('elbow: ', elbow[TIMES])) + 
	geom_errorbar(aes(x=numbers, ymin=wss-wss_std, ymax=wss+wss_std), width=0.25) +
	theme_minimal() +
	geom_vline(xintercept=elbow[TIMES]) 

#---- 11means ----

algo <- kmeans(PCA_data, centers=11)

aggregate(algo$cluster, by=list(algo$cluster), FUN=length)
data_kmeans <- list('data'=data, 'cluster'=algo$cluster)

fviz_cluster(data_kmeans, 
			 main = '11means',
			 frame.type = 'convex', 
			 geom = 'point', 
			 labelsize = 1,
) + theme_minimal() 

# ---- R squared ----

rsq <- function(data, cluster_label, r=2, target='data'){
	
	if(target=='cluster'){

		data <- as.matrix(data)
		n <- nrow(data)
		p <- ncol(data)
		
		
		for(i in 1:(ncol(data)-1)){
			RSS <- anova( aov( cluster_label ~ data ))[2,2] # RSS from linear regression
			TSS <- sum((cluster_label-mean(cluster_label))^2)
			rsq <- 1 - RSS/TSS * (n-1)/(n-p-1)
		}
		
		return(round(rsq, r))
	}
	
	
	if(target=='data'){
		
		data <- data.frame(data)
		data$cluster <- cluster_label
		rsq <- rep(NA, ncol(data)-1)
		
		for(i in 1:(ncol(data)-1)){
			RSS <- anova( aov( data[,i] ~ data$cluster))[2,2] # RSS from linear regression
			TSS <- sum((data[[i]])^2)
			rsq[i] <- 1 - RSS/TSS 
		}
		
		data <- data[, -ncol(data)]
		
	}
	
	rsq <- round(rsq, r)
	RSQ <- cbind(colnames(data), rsq) %>% data.frame(stringsAsFactors=FALSE)
	colnames(RSQ) <- c('variable', 'r_squared') 
	#RSQ <- RSQ[order(RSQ$r_squared), ]
	
	RSQ[dim(RSQ)[1]+1, ] <- c('mean', mean(as.numeric(RSQ[[2]])) %>% round(n)) # Add mean rsq 
	return(RSQ)
}

comp_data <- PCA_data

rsq_comparison <- merge(rsq(comp_data, sin_top), rsq(comp_data, cpl_top), 
						by.x='variable', by.y='variable') %>%
				  merge(rsq(comp_data, algo$cluster),
				  		by.x='variable', by.y='variable')
colnames(rsq_comparison) <- c('variable', 'rsq single linkage', 'rsq complete linkage', 'rsq 11means')
rsq_comparison

# ---- in R squared ----

comp_data <- PCA_data
rsq_inv_comp <- list('rsq single linkage' = rsq(comp_data, sin_top, target='cluster'), 
					 'rsq complete linkage' = rsq(comp_data, cpl_top, target='cluster'),
					 'rsq 11means' = rsq(comp_data, algo$cluster, target='cluster'))
rsq_inv_comp


comp_data <- data
rsq_inv_comp <- list('rsq single linkage' = rsq(comp_data, sin_top, target='cluster'), 
					 'rsq complete linkage' = rsq(comp_data, cpl_top, target='cluster'),
					 'rsq 11means' = rsq(comp_data, algo$cluster, target='cluster'))
rsq_inv_comp








