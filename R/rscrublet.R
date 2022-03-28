#' Example pbmc8k dataset
#'
#' A 8k pbmc single cell dataset
#'
#' @format dgTMatrix with UMI counts
#' @source \url{http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_filtered_gene_bc_matrices.tar.gz}
"pbmc8k"

#' Standard pipeline for preprocessing, doublet simulation, and doublet prediction
#' 
#' 
#' Automatically sets a threshold for calling doublets, but it's best to check 
#'        this by running \link{plot_doublet_histogram} afterwards and adjusting threshold manually
#'
#' @param E_obs sparse dgTMatrix (it will automatically coerced to dgTMatrix) n_cells*n_genes containing raw (unnormalized) UMI-based transcript counts. 
#' @param total_counts numerical vector of total UMI counts per cell. If NULL (default), this is calculated as the row sums of \code{E_obs}. 
#' @param sim_doublet_ratio Number of doublets to simulate relative to the number of observed transcriptomes.
#' @param n_neighbors Number of neighbors used to construct the KNN graph of observed transcriptomes and simulated doublets. If NULL (default), this is set to round(0.5 * sqrt(n_cells))
#' @param expected_doublet_rate The estimated doublet rate for the experiment.
#' @param stdev_doublet_rate Uncertainty in the expected doublet rate.
#' @param random_state Random seed for doublet simulation, approximate nearest neighbor search, and PCA/TruncatedSVD.
#' @param synthetic_doublet_umi_subsampling Rate for sampling UMIs when creating synthetic doublets. If 1.0, each doublet is created by simply adding the UMIs from two randomly sampled observed transcriptomes. For values less than 1, the UMI counts are added and then randomly sampled at the specified rate.
#' @param use_approx_neighbors logical (default TRUE). Use approximate nearest neighbor method (annoy) for the KNN classifier. Current implementation of exact KNN is much slower.
#' @param distance_metric Distance metric used when finding nearest neighbors. One of 'euclidean', 'manhattan', 'hamming', 'Angular' if use_approx_neighbors is true, or any method acceptable by \link{dist}
#' @param get_doublet_neighbor_parents logical, If TRUE, return the parent transcriptomes that generated the doublet neighbors of each observed transcriptome. This information can be used to infer the cell states that generated a given doublet state
#' @param min_counts  Used for gene filtering prior to PCA. Genes expressed at fewer than \code{min_counts} in fewer than \code{min_cells} (see below) are excluded.
#' @param min_cells Used for gene filtering prior to PCA. Genes expressed at fewer than \code{min_counts} (see above) in fewer than \code{min_cells} are excluded.
#' @param min_gene_variability_pctl Used for gene filtering prior to PCA. Keep the most highly variable genes (in the top min_gene_variability_pctl percentile), as measured by the v-statistic (Klein et al., Cell 2015).
#' @param log_transform logical (default: FALSE) If TRUE, log-transform the counts matrix (\code{log10(log_pseudocount+TPM)})
#' @param mean_center logical, If TURE (default), center the data such that each gene has a mean of 0.
#' @param normalize_variance logical, If TURE (default) If True, normalize the data such that each gene has a variance of 1.
#' @param n_prin_comps Number of principal components used to embed the transcriptomes prior to k-nearest-neighbor graph construction (default 30)
#' @param doublets_parents integer matrix with two columns that provide indexes for doublet parents. This option might be used to simulate doublets from specific parents, for instance for comparison with python version
#' @param log_pseudocount pseudocount for log transformation (default 1)
#' @param verbose If TRUE, print progress updates.
#'
#' @return list with following items:
#' \itemize{
#' \item doublet_scores_obs - Doublet scores for observed transcriptomes
#' \item doublet_scores_sim - Doublet scores for simulated doublets. 
#' \item doublet_errors_obs	- Standard error in the doublet scores for observed transcriptomes.
#' \item doublet_errors_sim	- Standard error in the doublet scores for simulated doublets.
#' \item doublet_neighbor_parent - parent transcriptomes that generated the doublet neighbors of each observed transcriptome. This information can be used to infer the cell states that generated a given doublet state
#' \item expected_doublet_rate - equal to \code{expected_doublet_rate} parameter
#' }
#' @export
#' 
#' @import RSpectra
#' @import RcppAnnoy
#' @import Matrix
#'
#' @examples
#' # run rscrublet of 8k pbmc example dataset
#' scrr = scrub_doublets(E_obs = pbmc8k,expected_doublet_rate=0.06,min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
#' # set threshould automatically 
#' scrr=call_doublets(scrr)
#' # examine score distribution
#' plot_doublet_histogram(scrr)
#' # find predicted doublets
#' rownames(pbmc8k)[scrr$predicted_doublets]
scrub_doublets = function(E_obs,
													total_counts=NULL, 
													sim_doublet_ratio=2.0, 
													n_neighbors=NULL, 
													expected_doublet_rate=0.1,
													stdev_doublet_rate=0.02,
													random_state=0,
													synthetic_doublet_umi_subsampling=1.0, 
													use_approx_neighbors=TRUE, 
													distance_metric='euclidean', 
													get_doublet_neighbor_parents=FALSE, 
													min_counts=3,
													min_cells=3, 
													min_gene_variability_pctl=85, 
													log_transform=FALSE, 
													mean_center=TRUE, 
													normalize_variance=TRUE,
													n_prin_comps=30,
													doublets_parents = NULL, # to use external doublets parent (to compare with python Scrablet for instance)
													log_pseudocount = 1,
													verbose=TRUE){
	# init
	if(verbose) cat('Preprocessing...','\n',sep='')
	time = Sys.time()
	E_obs = as(E_obs,'dgTMatrix')
	total_counts_obs = rowSums(E_obs)
	if(is.null(n_neighbors))
		n_neighbors = round(0.5*sqrt(nrow(E_obs))) 
	# norm
	if(verbose) cat(format(Sys.time()-time,digit=2),': Norm...','\n',sep='')
	E_obs_norm = tot_counts_norm(E_obs, total_counts=total_counts_obs)
	# filter
	if(verbose) cat(format(Sys.time()-time,digit=2),': Filter...','\n',sep='')
	gene_filter = filter_genes(E_obs_norm,
														 min_counts=min_counts,
														 min_cells=min_cells,
														 min_vscore_pctl=min_gene_variability_pctl)
	
	
	E_obs = E_obs[,gene_filter]
	E_obs_norm = E_obs_norm[,gene_filter]
	
	if(verbose) cat(format(Sys.time()-time,digit=2),': Simulating doublets...','\n',sep='')
	doub = simulate_doublets(E_obs,sim_doublet_ratio=sim_doublet_ratio,total_counts_obs=total_counts_obs,pair_ix = doublets_parents)
	
	E_sim = doub$E_sim
	total_counts_sim = doub$total_counts_sim
	
	E_obs_norm = tot_counts_norm(E_obs, total_counts=total_counts_obs,target_total = 1e6)
	E_sim_norm = tot_counts_norm(E_sim, total_counts=total_counts_sim,target_total = 1e6)
	
	
	
	if(log_transform){
		E_obs_norm = log_normalize(E_obs_norm, log_pseudocount)
		E_sim_norm = log_normalize(E_sim_norm, log_pseudocount)
	}
	
	if(normalize_variance){ 
		gene_stdevs = sqrt(sparse_col_var(E_obs_norm)) 
		E_obs_norm@x = E_obs_norm@x/gene_stdevs[E_obs_norm@j+1]
		E_sim_norm@x = E_sim_norm@x/gene_stdevs[E_sim_norm@j+1]
	}
	# matrices are still sparse 
	# and they will be sparse if mean_center is FALSE
	if(mean_center){
		gene_means = colMeans(E_obs_norm)
		E_obs_norm = sweep(as.array(E_obs_norm),2,gene_means,'-')
		E_sim_norm = sweep(as.array(E_sim_norm),2,gene_means,'-')
		
		# use PCA
		# pca = prcomp(E_obs_norm,center=FALSE,rank. = n_prin_comps) # consider to fins faster pca
		# manifold_obs = pca$x
		# manifold_sim = predict(pca,E_sim_norm)
		# RSpectra::svds seems to be much faster and theoretically (?) should be equal to PCA on centered data
		# Scrublet used truncated_svd in case of non-centered data (?) probably because in python PCA is faster than in R (?) and because svd (contrary to pca) doesnt work on sparse matrices
		# in R prcomp is very slow, so I switched to RSpectra::svds independently on mean_center. 
		# so in brief: original Scrublet used if(mean_center) pca else truncated_svd; I used svd in both cases
		
	}
	if(verbose) cat(format(Sys.time()-time,digit=2),': Embedding transcriptomes using SVD...','\n',sep='')
	svds = svds(E_obs_norm,k=n_prin_comps)
	manifold_obs = sweep(E_obs_norm %*% sweep(svds$v,2,svds$d,'/'),2,svds$d,'*')
	manifold_sim = sweep(E_sim_norm %*% sweep(svds$v,2,svds$d,'/'),2,svds$d,'*')
	
	if(verbose) cat(format(Sys.time()-time,digit=2),': Calculating doublet scores...','\n',sep='')
	nnc = nearest_neighbor_classifier(manifold_obs,
																		manifold_sim,
																		k=n_neighbors, 
																		use_approx_nn=use_approx_neighbors, 
																		distance_metric=distance_metric, 
																		exp_doub_rate=expected_doublet_rate, 
																		stdev_doub_rate=stdev_doublet_rate, 
																		get_neighbor_parents=get_doublet_neighbor_parents,
																		doub$doublet_parents,
																		random_state=random_state)
	if(verbose) cat('Elapsed time: ',format(Sys.time()-time,digit=2),'\n',sep='')
	nnc$expected_doublet_rate = expected_doublet_rate
	nnc
}

#' Automatically sets a threshold for calling doublets, 
#' but it's best to check this by running plot_doublet_histogram() afterwards and adjusting threshold by call_doublets with user defined \code{threshold} parameter
#'
#' @param nnc list, output of \link{scrub_doublets} function
#' @param threshold  Doublet score threshold for calling a transcriptome
#' a doublet. If NULL, this is set automatically by looking
#' for the minimum between the two modes of the \code{nnc$doublet_scores_sim}
#' histogram. It is best practice to check the threshold visually
#' using the \code{nnc$doublet_scores_sim} histogram and/or based on 
#' co-localization of predicted doublets in a 2-D embedding.
#' @param verbose  If TRUE, print summary statistics.
#'
#' @return \code{nnc} list with additional items:
#' #' \itemize{
#' \item predicted_doublets - logical index of predicted doublet cells
#' \item z_scores_ - Z-score conveying confidence in doublet calls. 
#' \item threshold - Doublet score threshold for calling a transcriptome a doublet.
#' \item detected_doublet_rate - Fraction of observed transcriptomes that have been called doublets.
#' \item detectable_doublet_fraction - Estimated fraction of doublets that are detectable, i.e., fraction of simulated doublets with doublet scores above \code{threshold}
#' \item overall_doublet_rate - Estimated overall doublet rate, \code{detected_doublet_rate/detectable_doublet_fraction}. Should agree (roughly) with \code{expected_doublet_rate}.
#' }
#' @export
#'
#' @examples
#' # run rscrublet of 8k pbmc example dataset
#' scrr = scrub_doublets(E_obs = pbmc8k,expected_doublet_rate=0.06,min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
#' # set threshould automatically 
#' scrr=call_doublets(scrr)
#' # examine score distribution
#' plot_doublet_histogram(scrr)
#' # find predicted doublets
#' rownames(pbmc8k)[scrr$predicted_doublets]
call_doublets = function(nnc,threshold=NULL,verbose=TRUE){
	if(is.null(threshold)){ 
		d = density(nnc$doublet_scores_sim,from=min(nnc$doublet_scores_sim),to=max(nnc$doublet_scores_sim),bw = 0.05)
		l = length(d$y) 
		maxs = which(d$y>c(0,d$y[-l]) & d$y>c(d$y[-1],0))
		mins = which(d$y<c(0,d$y[-l]) & d$y<c(d$y[-1],0))
		if(length(maxs)>1 & length(maxs)>0)
			mins = mins[mins>maxs[1] & mins <maxs[2]]
		
		if(length(mins) != 1){
			message('Fail to set threshold automatically, 0.2 is used by default. Rerun `call_doublets` with user-specified threshold.')
			threshold = 0.2
		}else{
			threshold = d$x[mins]
			if(verbose) cat('Automatically set threshold at doublet score ',format(threshold,digits=4),'\n',sep='')
		}
	}
	
	
	Ld_obs = nnc$doublet_scores_obs
	Ld_sim = nnc$doublet_scores_sim
	se_obs = nnc$doublet_errors_obs
	Z = (Ld_obs - threshold) / se_obs
	nnc$predicted_doublets = Ld_obs > threshold
	nnc$z_scores_ = Z
	nnc$threshold = threshold
	nnc$detected_doublet_rate = sum(Ld_obs>threshold) / length(Ld_obs)
	nnc$detectable_doublet_fraction = sum(Ld_sim>threshold) / length(Ld_sim)
	nnc$overall_doublet_rate = nnc$detected_doublet_rate / nnc$detectable_doublet_fraction
	
	if(verbose){
		cat('Detected doublet rate = ',format(100*nnc$detected_doublet_rate,digits=2),'%\n')
		cat('Estimated detectable doublet fraction = ',format(100*nnc$detectable_doublet_fraction,digits=2),'%\n')
		cat('Overall doublet rate:\n')
		cat('\tExpected   = ',format(100*nnc$expected_doublet_rate,digits=2),'%\n')
		cat('\tEstimated  = ',format(100*nnc$overall_doublet_rate,digits=2),'%\n')
	}
	nnc
}

#' Plot histogram of doublet scores for observed transcriptomes and simulated doublets 
#'
#' The histogram for simulated doublets is useful for determining the correct doublet 
#' score threshold. To set threshold to a specific value, T, run call_doublets(threshold=thr.value).
#'
#' @param nnc list, output of call_doublets function
#' @param breaks number of breaks to be used in hist
#'
#' @export
#'
#' @examples
#' #' scrr = scrub_doublets(E_obs = pbmc8k,expected_doublet_rate=0.06,min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
#' scrr=call_doublets(scrr)
#' plot_doublet_histogram(scrr)
plot_doublet_histogram = function(nnc,breaks=30){
	par(mfrow=c(1,2),tck=-0.02,mgp=c(1.1,0.2,0),mar=c(3,3,1.5,0),oma=c(0,0,0,1))
	h=hist(nnc$doublet_scores_obs,breaks =breaks,plot=F)
	h$counts = log(h$counts+1)
	plot(h,xlab='Doublet score',ylab='Prob. density',main='Observed Transcriptomes',border = NA,col='gray',yaxt='n')
	abline(v=nnc$threshold,col='red')
	lab = 10^(0:10)
	at = log(lab+1)
	lab = lab[at < max(h$counts)]
	at = at[at < max(h$counts)]
	axis(2,at,lab)
	hist(nnc$doublet_scores_sim,breaks=breaks,xlab='Doublet score',ylab='Prob. density',main='Simulated doublets',border = NA,col='gray')
	abline(v=nnc$threshold,col='red')
}


tot_counts_norm = function(E, total_counts = NULL, exclude_dominant_frac = 1, included = NULL, target_total = NULL){
	# Cell-level total counts normalization of input counts matrix, excluding overly abundant genes if desired.
	# Return normalized counts, average total counts, and (if exclude_dominant_frac < 1) list of genes used to calculate total counts 
	if(is.null(total_counts)){
		if(is.null(included) == 0){
			if(exclude_dominant_frac == 1){
				tots_use = rowSums(E)
			}else{
				En = E
				En@x = En@x/rowSums(En)[En@i+1]
				included = colSums(En>exclude_dominant_frac)==0
				tots_use = rowSums(E[,included])
				print(paste0('Excluded ',sum(!included),' genes from normalization'))
			}
		}else{
			tots_use = rowSums(E[,included])
		}
	}else{
		tots_use = total_counts
	}
	
	if(is.null(target_total))
		target_total = mean(tots_use)
	En = E
	En@x = En@x/tots_use[En@i+1]*target_total
	En
}

runningquantile = function(x, y, p, nBins){
	ind = order(x)
	x = x[ind]
	y = y[ind]
	
	dx = (x[length(x)] - x[1]) / nBins
	
	xOut = seq(x[1]+dx/2, x[length(x)]-dx/2, length.out= nBins)
	
	yOut = rep(0,length(xOut))
	
	for(i in 1:length(xOut)){
		ind = which((x >= xOut[i]-dx/2) & (x < xOut[i]+dx/2)) # can be optimized
		if(length(ind) > 0){
			yOut[i] = quantile(y[ind], p/100) # so small p looks weird, but that is how it work in Scrublet. it might be a mistake (0.1% instead of probability=0.1), btw i would say that median looks more reasonable
		}else{
			if(i > 1)
				yOut[i] = yOut[i-1]
			else
				yOut[i] = NaN
		}
	}
	list(x=xOut, y=yOut)
}

get_vscores = function(E, min_mean=0, nBins=50, fit_percentile=0.1, error_wt=1){
	
	## Calculate v-score (above-Poisson noise statistic) for genes in the input counts matrix
	## Return v-scores and other stats
	
	ncell = nrow(E)
	
	mu_gene = colMeans(E)
	gene_ix = which(mu_gene > min_mean)
	mu_gene = mu_gene[gene_ix]
	
	tmp = E[,gene_ix]
	
	tmp@x = tmp@x^2
	
	var_gene = colMeans(tmp) - mu_gene ^ 2
	
	rm(tmp)
	FF_gene = var_gene / mu_gene
	
	data_x = log(mu_gene)
	data_y = log(FF_gene / mu_gene)
	
	xy = runningquantile(data_x, data_y, fit_percentile, nBins)
	
	x = xy$x[!is.na(xy$y)]
	y = xy$y[!is.na(xy$y)]
	rm(xy)
	
	t = log(FF_gene[mu_gene>0])
	b = seq(min(t),max(t),length.out=201)
	h = hist(t, breaks = b,plot = F)$counts
	
	b = b[-length(b)] + diff(b)/2
	max_ix = which.max(h)
	c = max(exp(b[max_ix]), 1)
	
	gLog = function(input) log(input[[2]] * exp(-input[[1]]) + input[[3]])
	errFun = function(b2)sum(abs(gLog(list(x,c,b2))-y) ** error_wt)
	b = optimize(f = errFun, interval = c(0,100))$minimum # results are slightly different from python, likely due to tol
	a = c / (1+b) - 1
	
	
	v_scores = FF_gene / ((1+a)*(1+b) + b * mu_gene);
	CV_eff = sqrt((1+a)*(1+b) - 1);
	CV_input = sqrt(b);
	
	gc(verbose = FALSE)
	list(v_scores=v_scores, 
			 CV_eff=CV_eff,
			 CV_input=CV_input,
			 gene_ix=gene_ix,
			 mu_gene=mu_gene,
			 FF_gene=FF_gene, 
			 a=a,
			 b=b)
}


filter_genes = function(E, base_ix = NULL, min_vscore_pctl = 85, min_counts = 3, min_cells = 3, show_vscore_plot = FALSE, sample_name = ''){
	## Filter genes by expression level and variability
	## Return list of filtered gene indices
	
	if(is.null(base_ix))
		base_ix = 1:nrow(E)
	
	scores = get_vscores(E[base_ix, ])
	ix2 = scores$v_scores>0
	Vscores = scores$v_scores[ix2]
	gene_ix = scores$gene_ix[ix2]
	mu_gene = scores$mu_gene[ix2]
	FF_gene = scores$FF_gene[ix2]
	min_vscore = quantile(Vscores, min_vscore_pctl/100)
	
	ix = ((colSums(E[,gene_ix] >= min_counts) >= min_cells) & (Vscores >= min_vscore))
	
	if(show_vscore_plot){
		message("show_vscore_plot is not implemented yet")
	}
	gene_ix[ix]
}

subsample_counts = function(E, rate, original_totals, random_seed=0){
	if(rate < 1){
		set.seed(random_seed)
		E@x = rbinom(length(E@x),E@x, rate)
		current_totals = rowSums(E)
		unsampled_orig_totals = original_totals - current_totals
		unsampled_downsamp_totals = rbinom(length(unsampled_orig_totals),unsampled_orig_totals, rate)
		final_downsamp_totals = current_totals + unsampled_downsamp_totals
	}else
		final_downsamp_totals = original_totals
	list(E_sim=E, total_counts_sim=final_downsamp_totals)
}

simulate_doublets = function(E_obs,sim_doublet_ratio, synthetic_doublet_umi_subsampling=1.0,random_state=0,total_counts_obs,pair_ix=NULL){
	n_obs = nrow(E_obs)
	n_sim = as.integer(n_obs * sim_doublet_ratio)
	
	set.seed(random_state)
	if(is.null(pair_ix))
		pair_ix = matrix(sample(n_obs,n_sim*2,replace = TRUE),ncol=2)
	
	E1 = E_obs[pair_ix[,1],]
	E2 = E_obs[pair_ix[,2],]
	tots1 = total_counts_obs[pair_ix[,1]]
	tots2 = total_counts_obs[pair_ix[,2]]
	if(synthetic_doublet_umi_subsampling < 1){
		res = subsample_counts(E1+E2, synthetic_doublet_umi_subsampling, tots1+tots2, random_seed=random_state)
	}else
		res = list(E_sim = E1+E2,total_counts_sim = tots1+tots2)
	res$doublet_parents = pair_ix
	res
}

log_normalize = function(X,pseudocount=1){
	X@x = log10(X@x + pseudocount)
	X
}

sparse_col_var = function(x){
	v = split(x@x,x@j)
	v = sapply(v,function(z){
		m = sum(z)/nrow(x)
		(sum((z - m)^2)+(m^2)*(nrow(x)-length(z)))/(nrow(x)-1)
	})
	r = setNames(rep(NA_real_,ncol(x)),colnames(x))
	r[colnames(x)[as.numeric(names(v))+1]] = v
	r
}

get_knn_graph = function(X, k=5, dist_metric='euclidean', approx=FALSE, return_edges=TRUE, random_seed=0){
	if(approx){
		a = new(get(paste0('Annoy',first2Upper(dist_metric))), ncol(X))
		a$setSeed(random_seed)
		for(i in 1:nrow(X))
			a$addItem(i, X[i,])
		
		a$build(10) # 10 trees
		
		knn = matrix(NA,nrow=nrow(X),ncol=k)
		for(iCell in 1:nrow(X))
			knn[iCell,] = a$getNNsByItem(iCell,k)
		
	}else{
		d = as.matrix(dist(X,method = dist_metric))
		knn = apply(d,1,function(x)order(x)[1:(k+1)])
		knn = t(knn)
		knn = t(sapply(1:nrow(knn),function(i){if(knn[i,1]==i)knn[i,-1]else setdiff(knn[i,],i)}))
		rm(d);gc()
	}
	knn
}

nearest_neighbor_classifier = function(manifold_obs, manifold_sim,k=40, use_approx_nn=TRUE, distance_metric='euclidean', 
																			 exp_doub_rate=0.1, stdev_doub_rate=0.03, get_neighbor_parents=FALSE,doublet_parents,random_state){
	manifold = rbind(manifold_obs, manifold_sim)
	
	n_obs = nrow(manifold_obs)
	n_sim = nrow(manifold_sim)
	
	doub_labels = c(rep(0,n_obs), rep(1,n_sim))
	
	# Adjust k (number of nearest neighbors) based on the ratio of simulated to observed cells
	k_adj = as.integer(round(k * (1+n_sim/n_obs)))

	# Find k_adj nearest neighbors
	neighbors = get_knn_graph(manifold, k=k_adj,dist_metric=distance_metric, approx=use_approx_nn, return_edges=FALSE, random_seed=random_state)
	
	# Calculate doublet score based on ratio of simulated cell neighbors vs. observed cell neighbors
	n_sim_neigh = apply(neighbors, 1, function(x)sum(doub_labels[x]))
	n_obs_neigh = ncol(neighbors) - n_sim_neigh
	
	rho = exp_doub_rate
	r = n_sim / n_obs
	nd = n_sim_neigh
	ns = n_obs_neigh
	N = k_adj
	
	# Bayesian
	q=(nd+1)/(N+2)
	Ld = q*rho/r/(1-rho-q*(1-rho-rho/r))
	
	se_q = sqrt(q*(1-q)/(N+3))
	se_rho = stdev_doub_rate
	
	se_Ld = q*rho/r / (1-rho-q*(1-rho-rho/r))**2 * sqrt((se_q/q*(1-rho))**2 + (se_rho/rho*(1-q))**2)
	
	doublet_scores_obs = Ld[doub_labels == 0]
	doublet_scores_sim = Ld[doub_labels == 1]
	doublet_errors_obs = se_Ld[doub_labels==0]
	doublet_errors_sim = se_Ld[doub_labels==1]
	
	# get parents of doublet neighbors, if requested
	neighbor_parents = NULL
	if(get_neighbor_parents){
		parent_cells = doublet_parents
		neighbors = neighbors - n_obs
		neighbor_parents = list()
		for(iCell in 1:n_obs){
			this_doub_neigh = neighbors[iCell,][neighbors[iCell,] > 0]
			if(length(this_doub_neigh) > 0){
				this_doub_neigh_parents = unique(as.integer(parent_cells[this_doub_neigh,]))
				neighbor_parents[[length(neighbor_parents)+1]] = this_doub_neigh_parents
			}else
				neighbor_parents[length(neighbor_parents)+1] = list(NULL)
		}
	}
	list(doublet_scores_obs=doublet_scores_obs,
			 doublet_scores_sim=doublet_scores_sim,
			 doublet_errors_obs=doublet_errors_obs,
			 doublet_errors_sim=doublet_errors_sim,
			 doublet_neighbor_parent=neighbor_parents)
}

first2Upper = function(t){
	paste0(toupper(substr(t,1,1)),substr(t,2,nchar(t)))
}