# Center
center_scale = function(x) {
  scale(x, scale = FALSE)
}

# Plot function
ggplot.prep = function(pvalue){
  ecount = c()
  pcount = c()
  n = length(pvalue)
  for (i in 1:n) {
    critical.value = i/n
    ind = which(pvalue<critical.value)
    ecount[i] = length(ind)
  }
  x = seq(0,1,1/n)
  x = x[-1]
  ecount = ecount/n
  type.I.error = data.frame(Critical.value = x, P.value=ecount)
  return(type.I.error)
}

SKAT.pseudobulk = function(geno, annotation) {
  counts.mat = geno[[1]]
  cell.types = levels(annotation$cell.type)
  donors = levels(annotation$donor)
  n.donors = length(donors)
  n.genes = nrow(counts.mat)
  n.cell.types = length(cell.types)

  # Initialize list to store pseudobulk matrices
  pseudobulk.list = list()
  # Process each cell type
  for (ct in cell.types) {
    # Initialize pseudobulk matrix for this cell type (n.donors x n.genes)
    pseudobulk.mat = matrix(0,nrow = n.donors, ncol = n.genes,dimnames = list(donors, rownames(counts.mat)))

    # Process each donor
    for (i in 1:n.donors) {
      donor = donors[i]
      # Get indices of cells for this donor and cell type
      cell.idx = which(annotation$donor == donor & annotation$cell.type == ct)

      if (length(cell.idx) > 0) {
        # Calculate mean expression across cells
        # If only one cell, just take that column; otherwise, take rowMeans
        if (length(cell.idx) == 1) {
           pseudobulk.mat[i, ] = counts.mat[, cell.idx]
        } else {
            pseudobulk.mat[i, ] = rowMeans(counts.mat[, cell.idx, drop = FALSE])
        }
      }
    }
    pseudobulk.mat = log2(pseudobulk.mat + 1)
    pseudobulk.mat = scale(pseudobulk.mat, center = TRUE, scale = TRUE)
    # Handle genes with zero variance (set to 0)
    zero.var = which(is.na(pseudobulk.mat[1, ]))
    if (length(zero.var) > 0) pseudobulk.mat[, zero.var] = 0
    # Add to list with cell type name
    pseudobulk.list[[ct]] = pseudobulk.mat
  }
  X.SKAT = matrix(unlist(pseudobulk.list),n.donors,(n.cell.types*n.genes))
  return(X.SKAT)

}

# Source package
library(glmnet)
library(MASS)
library(matrixcalc)
library(bracer)
library(CompQuadForm)
library(expm)
library(ggpubr)
library(ggplot2)
library(purrr)
library(optparse)
library(Seurat)
library(nebula)
library(MAST)
library(SKAT)

option.list = list(
  make_option(c("-s", "--simfun"), type="character", default="linear",
              help="Simulation nonlinear function", metavar="character"),
  make_option(c("-p", "--para"), type="character", default='500_4_1',
              help="SNP and scale parameter", metavar="character"),
  make_option(c("-t", "--times"), type="integer", default=200,
              help="Simulation times", metavar="character"),
  make_option(c("-l", "--lambda"), type="double", default=0,
              help="Penalized parameter", metavar="character"))
opt.parser = OptionParser(option_list=option.list)
opt = parse_args(opt.parser)
simfun = opt$simfun
para = opt$para
max.time = opt$times
lambda = opt$lambda
#####################
#simfun = 'linear'
para = '10_0_0.2'
max.time = 200
lambda = 0.5
#####################

times.list = 1:max.time
#lambda.list = c(0,0.05,0.5,10,50)
first.order.ind = 3
second.order.ind= 6

# Source functions
path = "/home/heng.ge/blue/heng.ge/KNN_TWAS_Test_Simulation/"
function.path = paste0(path,'R/')
data.path = paste0(path,'data/')
save.path = paste0(path,'result/')
walk(list.files(function.path,pattern="*.R$",full.names = T),source,.GlobalEnv)

# Initialize p-value vectors for PKNN
pvalue.PKNN = c()
pvalue.PKNN.first = c()
pvalue.PKNN.second = c()

pvalue.PKNN.null = c()
pvalue.PKNN.first.null = c()
pvalue.PKNN.second.null = c()

# Initialize p-value vectors for NEBULA
pvalue.NEBULA = c()
pvalue.NEBULA.null = c()

# Initialize p-value vectors for MAST
pvalue.MAST = c()
pvalue.MAST.null = c()

pvalue.SKAT = c()
pvalue.SKAT.null = c()

for (times in times.list){

  # Set path
  sim.path = paste0(data.path,simfun,'/')
  sim.name = paste0('Sim_',times,'.Rdata')
  setwd(sim.path)
  load(paste0(sim.path,sim.name))
  X.list = list(counts.mat)

  # Read data and split data
  y = pheno.donor
  y.null = rnorm(length(y))

  #================== PKNN Analysis ==================
  # Estimate parameter
  input.kernel = rep(c('product'),3)


  input.kernel = rep(list(c('product')),3)
  result.PKNN = PKNN.mixchiq.overall.test(y,X.list,test.ind=c(first.order.ind,second.order.ind),input.kernel=input.kernel,lambda = lambda,annotation=annotation,KNN.type = 'cell')
  result.PKNN.first = PKNN.mixchiq.ind.test(y,X.list,first.order.ind,input.kernel,lambda=lambda,annotation=annotation,KNN.type = 'cell')
  result.PKNN.second = PKNN.mixchiq.ind.test(y,X.list,second.order.ind,input.kernel,lambda=lambda,annotation=annotation,KNN.type = 'cell')

  result.PKNN.null = PKNN.mixchiq.overall.test(y.null,X.list,test.ind=c(first.order.ind,second.order.ind),input.kernel=input.kernel,lambda=lambda,annotation=annotation,KNN.type = 'cell')
  result.PKNN.first.null = PKNN.mixchiq.ind.test(y.null,X.list,first.order.ind,input.kernel,lambda=lambda,annotation=annotation,KNN.type = 'cell')
  result.PKNN.second.null = PKNN.mixchiq.ind.test(y.null,X.list,second.order.ind,input.kernel,lambda=lambda,annotation=annotation,KNN.type = 'cell')

  # Results for PKNN
  pvalue.PKNN = c(pvalue.PKNN,result.PKNN$pvalue)
  pvalue.PKNN.first = c(pvalue.PKNN.first,result.PKNN.first$pvalue)
  pvalue.PKNN.second = c(pvalue.PKNN.second,result.PKNN.second$pvalue)

  pvalue.PKNN.null = c(pvalue.PKNN.null,result.PKNN.null$pvalue)
  pvalue.PKNN.first.null = c(pvalue.PKNN.first.null,result.PKNN.first.null$pvalue)
  pvalue.PKNN.second.null = c(pvalue.PKNN.second.null,result.PKNN.second.null$pvalue)

  #================== SKAT Analysis ==================
  X.SKAT = SKAT.pseudobulk(X.list, annotation)

  obj = SKAT_Null_Model(y ~ 1, out_type = "C")
  obj.null = SKAT_Null_Model(y.null ~ 1, out_type = "C")

  result.SKAT = SKAT(Z = X.SKAT, obj = obj)
  result.SKAT.null = SKAT(Z = X.SKAT, obj = obj.null)

  pvalue.SKAT = c(pvalue.SKAT,result.SKAT$p.value)
  pvalue.SKAT.null = c(pvalue.SKAT.null,result.SKAT.null$p.value)

  #================== NEBULA Analysis ==================
  # Create design matrix with continuous phenotype
  pred_data = data.frame(
    cont_pheno = y[annotation$donor],
    row.names = colnames(counts.mat)
  )
  pred_nebula = model.matrix(~ cont_pheno, data = pred_data)
  # Fit NEBULA
  res_nebula = nebula(
    count  = counts.mat,
    id     = annotation$donor,  # Donor ID for each cell
    pred   = pred_nebula,
    offset = log(colSums(counts.mat)),
    method = "LN"
  )
  nebula_pvals = res_nebula$summary[, "p_cont_pheno"]
  nebula_pvals_clean = nebula_pvals[!is.na(nebula_pvals) & nebula_pvals > 0 & nebula_pvals <= 1]
  fisher_stat = -2 * sum(log(nebula_pvals_clean))
  df_fisher = 2 * length(nebula_pvals_clean)
  combined_pval = pchisq(fisher_stat, df = df_fisher, lower.tail = FALSE)
  pvalue.NEBULA = c(pvalue.NEBULA, combined_pval)

  #================== NEBULA NULL Analysis ==================
  pred_data = data.frame(
    cont_pheno = y.null[annotation$donor],
    row.names = colnames(counts.mat)
  )
  pred_nebula = model.matrix(~ cont_pheno, data = pred_data)
  # Fit NEBULA
  res_nebula = nebula(
    count  = counts.mat,
    id     = annotation$donor,  # Donor ID for each cell
    pred   = pred_nebula,
    offset = log(colSums(counts.mat)),
    method = "LN"
  )
  nebula_pvals = res_nebula$summary[, "p_cont_pheno"]
  nebula_pvals_clean = nebula_pvals[!is.na(nebula_pvals) & nebula_pvals > 0 & nebula_pvals <= 1]
  fisher_stat = -2 * sum(log(nebula_pvals_clean))
  df_fisher = 2 * length(nebula_pvals_clean)
  combined_pval = pchisq(fisher_stat, df = df_fisher, lower.tail = FALSE)
  pvalue.NEBULA.null = c(pvalue.NEBULA.null, combined_pval)

  #================== MAST Analysis ==================
  # Prepare data for MAST
  logcounts = log2(counts.mat + 1)
  ncells = colSums(counts.mat > 0)
  meta = data.frame(cont_pheno = y[annotation$donor],ncells = ncells,sid = as.numeric(annotation$donor),row.names = colnames(counts.mat))
  fdata = data.frame(primerid = rownames(counts.mat),row.names = rownames(counts.mat),stringsAsFactors = FALSE)
  sca = FromMatrix(exprsArray = as.matrix(logcounts),cData = meta,fData = fdata)
  zlm_mod = zlm(formula = ~ cont_pheno + ncells,sca = sca)
  # Likelihood ratio test
  lrt_result = lrTest(zlm_mod, "cont_pheno")
  # Extract hurdle p-values
  p_values_hurdle = lrt_result[, "hurdle", "Pr(>Chisq)"]
  p_values_hurdle_clean = p_values_hurdle[!is.na(p_values_hurdle) & p_values_hurdle > 0 & p_values_hurdle <= 1]
  fisher_stat_mast = -2 * sum(log(p_values_hurdle_clean))
  df_fisher_mast = 2 * length(p_values_hurdle_clean)
  combined_pval_mast = pchisq(fisher_stat_mast, df = df_fisher_mast, lower.tail = FALSE)
  pvalue.MAST = c(pvalue.MAST, combined_pval_mast)

  #================== MAST NULL Analysis ==================
  # Prepare data for MAST
  logcounts = log2(counts.mat + 1)
  ncells = colSums(counts.mat > 0)
  meta = data.frame(cont_pheno = y.null[annotation$donor],ncells = ncells,sid = as.numeric(annotation$donor),row.names = colnames(counts.mat))
  fdata = data.frame(primerid = rownames(counts.mat),row.names = rownames(counts.mat),stringsAsFactors = FALSE)
  sca = FromMatrix(exprsArray = as.matrix(logcounts),cData = meta,fData = fdata)
  zlm_mod = zlm(formula = ~ cont_pheno + ncells,sca = sca)
  # Likelihood ratio test
  lrt_result = lrTest(zlm_mod, "cont_pheno")
  # Extract hurdle p-values
  p_values_hurdle = lrt_result[, "hurdle", "Pr(>Chisq)"]
  p_values_hurdle_clean = p_values_hurdle[!is.na(p_values_hurdle) & p_values_hurdle > 0 & p_values_hurdle <= 1]
  fisher_stat_mast = -2 * sum(log(p_values_hurdle_clean))
  df_fisher_mast = 2 * length(p_values_hurdle_clean)
  combined_pval_mast = pchisq(fisher_stat_mast, df = df_fisher_mast, lower.tail = FALSE)
  pvalue.MAST.null = c(pvalue.MAST.null, combined_pval_mast)

  print('==========================Power==========================')
  print('==========================PKNN=========================')
  print(result.PKNN)
  print(result.PKNN.first)
  print(result.PKNN.second)
  print('==========================NEBULA=======================')
  print(paste("NEBULA p-value:", tail(pvalue.NEBULA, 1)))
  print('==========================MAST=========================')
  print(paste("MAST p-value:", tail(pvalue.MAST, 1)))
  print('=======================Type I error======================')
  print('==========================PKNN=========================')
  print(result.PKNN.null)
  print(result.PKNN.first.null)
  print(result.PKNN.second.null)
  print('==========================NEBULA=======================')
  print(paste("NEBULA null p-value:", tail(pvalue.NEBULA.null, 1)))
  print('==========================MAST=========================')
  print(paste("MAST null p-value:", tail(pvalue.MAST.null, 1)))

  save(pvalue.PKNN,pvalue.PKNN.first,pvalue.PKNN.second,
       pvalue.PKNN.null,pvalue.PKNN.first.null,pvalue.PKNN.second.null,
       pvalue.NEBULA, pvalue.NEBULA.null,
       pvalue.MAST, pvalue.MAST.null,pvalue.SKAT,pvalue.SKAT.null,
       file=paste0(save.path,simfun,'_',para,'_mixchiq_with_comparisons.Rdata'))

}

# Create plots for all methods
# PKNN plots
type.I.error = ggplot.prep(pvalue.SKAT)
p1 = ggplot(type.I.error, aes(x=Critical.value, y=P.value)) + xlim(0,1) + ylim(0,1)+
  geom_point(size=1)+geom_text(aes(label=ifelse(Critical.value==0.05,as.character(P.value),'')),hjust=0,vjust=1.5) +
  geom_abline(intercept = 0, slope = 1,color='red')+geom_vline(xintercept = 0.05, linetype="dotted",color = "blue")+
  geom_hline(yintercept = 0.05, linetype="dotted",color = "blue") + ggtitle('pvalue.SKAT') +
  theme(plot.title = element_text(hjust = 0.5))

type.I.error = ggplot.prep(pvalue.SKAT.null)
p2 = ggplot(type.I.error, aes(x=Critical.value, y=P.value)) + xlim(0,1) + ylim(0,1)+
  geom_point(size=1)+geom_text(aes(label=ifelse(Critical.value==0.05,as.character(P.value),'')),hjust=0,vjust=1.5) +
  geom_abline(intercept = 0, slope = 1,color='red')+geom_vline(xintercept = 0.05, linetype="dotted",color = "blue")+
  geom_hline(yintercept = 0.05, linetype="dotted",color = "blue") + ggtitle('pvalue.SKAT.null') +
  theme(plot.title = element_text(hjust = 0.5))

type.I.error = ggplot.prep(pvalue.PKNN)
p3 = ggplot(type.I.error, aes(x=Critical.value, y=P.value)) + xlim(0,1) + ylim(0,1)+
  geom_point(size=1)+geom_text(aes(label=ifelse(Critical.value==0.05,as.character(P.value),'')),hjust=0,vjust=1.5) +
  geom_abline(intercept = 0, slope = 1,color='red')+geom_vline(xintercept = 0.05, linetype="dotted",color = "blue")+
  geom_hline(yintercept = 0.05, linetype="dotted",color = "blue") + ggtitle('pvalue.PKNN') +
  theme(plot.title = element_text(hjust = 0.5))

type.I.error = ggplot.prep(pvalue.PKNN.first)
p6 = ggplot(type.I.error, aes(x=Critical.value, y=P.value)) + xlim(0,1) + ylim(0,1)+
  geom_point(size=1)+geom_text(aes(label=ifelse(Critical.value==0.05,as.character(P.value),'')),hjust=0,vjust=1.5) +
  geom_abline(intercept = 0, slope = 1,color='red')+geom_vline(xintercept = 0.05, linetype="dotted",color = "blue")+
  geom_hline(yintercept = 0.05, linetype="dotted",color = "blue") + ggtitle('pvalue.PKNN.first') +
  theme(plot.title = element_text(hjust = 0.5))

type.I.error = ggplot.prep(pvalue.PKNN.second)
p7 = ggplot(type.I.error, aes(x=Critical.value, y=P.value)) + xlim(0,1) + ylim(0,1)+
  geom_point(size=1)+geom_text(aes(label=ifelse(Critical.value==0.05,as.character(P.value),'')),hjust=0,vjust=1.5) +
  geom_abline(intercept = 0, slope = 1,color='red')+geom_vline(xintercept = 0.05, linetype="dotted",color = "blue")+
  geom_hline(yintercept = 0.05, linetype="dotted",color = "blue") + ggtitle('pvalue.PKNN.second') +
  theme(plot.title = element_text(hjust = 0.5))

type.I.error = ggplot.prep(pvalue.PKNN.null)
p10 = ggplot(type.I.error, aes(x=Critical.value, y=P.value)) + xlim(0,1) + ylim(0,1)+
  geom_point(size=1)+geom_text(aes(label=ifelse(Critical.value==0.05,as.character(P.value),'')),hjust=0,vjust=1.5) +
  geom_abline(intercept = 0, slope = 1,color='red')+geom_vline(xintercept = 0.05, linetype="dotted",color = "blue")+
  geom_hline(yintercept = 0.05, linetype="dotted",color = "blue") + ggtitle('pvalue.PKNN.null') +
  theme(plot.title = element_text(hjust = 0.5))

type.I.error = ggplot.prep(pvalue.PKNN.first.null)
p13 = ggplot(type.I.error, aes(x=Critical.value, y=P.value)) + xlim(0,1) + ylim(0,1)+
  geom_point(size=1)+geom_text(aes(label=ifelse(Critical.value==0.05,as.character(P.value),'')),hjust=0,vjust=1.5) +
  geom_abline(intercept = 0, slope = 1,color='red')+geom_vline(xintercept = 0.05, linetype="dotted",color = "blue")+
  geom_hline(yintercept = 0.05, linetype="dotted",color = "blue") + ggtitle('pvalue.PKNN.first.null') +
  theme(plot.title = element_text(hjust = 0.5))

type.I.error = ggplot.prep(pvalue.PKNN.second.null)
p14 = ggplot(type.I.error, aes(x=Critical.value, y=P.value)) + xlim(0,1) + ylim(0,1)+
  geom_point(size=1)+geom_text(aes(label=ifelse(Critical.value==0.05,as.character(P.value),'')),hjust=0,vjust=1.5) +
  geom_abline(intercept = 0, slope = 1,color='red')+geom_vline(xintercept = 0.05, linetype="dotted",color = "blue")+
  geom_hline(yintercept = 0.05, linetype="dotted",color = "blue") + ggtitle('pvalue.PKNN.second.null') +
  theme(plot.title = element_text(hjust = 0.5))

# NEBULA plots
type.I.error = ggplot.prep(pvalue.NEBULA)
p_nebula = ggplot(type.I.error, aes(x=Critical.value, y=P.value)) + xlim(0,1) + ylim(0,1)+
  geom_point(size=1)+geom_text(aes(label=ifelse(Critical.value==0.05,as.character(P.value),'')),hjust=0,vjust=1.5) +
  geom_abline(intercept = 0, slope = 1,color='red')+geom_vline(xintercept = 0.05, linetype="dotted",color = "blue")+
  geom_hline(yintercept = 0.05, linetype="dotted",color = "blue") + ggtitle('pvalue.NEBULA') +
  theme(plot.title = element_text(hjust = 0.5))

type.I.error = ggplot.prep(pvalue.NEBULA.null)
p_nebula_null = ggplot(type.I.error, aes(x=Critical.value, y=P.value)) + xlim(0,1) + ylim(0,1)+
  geom_point(size=1)+geom_text(aes(label=ifelse(Critical.value==0.05,as.character(P.value),'')),hjust=0,vjust=1.5) +
  geom_abline(intercept = 0, slope = 1,color='red')+geom_vline(xintercept = 0.05, linetype="dotted",color = "blue")+
  geom_hline(yintercept = 0.05, linetype="dotted",color = "blue") + ggtitle('pvalue.NEBULA.null') +
  theme(plot.title = element_text(hjust = 0.5))

# MAST plots
type.I.error = ggplot.prep(pvalue.MAST)
p_mast = ggplot(type.I.error, aes(x=Critical.value, y=P.value)) + xlim(0,1) + ylim(0,1)+
  geom_point(size=1)+geom_text(aes(label=ifelse(Critical.value==0.05,as.character(P.value),'')),hjust=0,vjust=1.5) +
  geom_abline(intercept = 0, slope = 1,color='red')+geom_vline(xintercept = 0.05, linetype="dotted",color = "blue")+
  geom_hline(yintercept = 0.05, linetype="dotted",color = "blue") + ggtitle('pvalue.MAST') +
  theme(plot.title = element_text(hjust = 0.5))

type.I.error = ggplot.prep(pvalue.MAST.null)
p_mast_null = ggplot(type.I.error, aes(x=Critical.value, y=P.value)) + xlim(0,1) + ylim(0,1)+
  geom_point(size=1)+geom_text(aes(label=ifelse(Critical.value==0.05,as.character(P.value),'')),hjust=0,vjust=1.5) +
  geom_abline(intercept = 0, slope = 1,color='red')+geom_vline(xintercept = 0.05, linetype="dotted",color = "blue")+
  geom_hline(yintercept = 0.05, linetype="dotted",color = "blue") + ggtitle('pvalue.MAST.null') +
  theme(plot.title = element_text(hjust = 0.5))

# Combined plot with all methods
pp_all = ggarrange(p3,p6,p7,p1,p_nebula,p_mast,p10,p13,p14,p2,p_nebula_null,p_mast_null,nrow=2,ncol=6)
png(paste0(save.path,simfun,'_',para,'.png'),width = 1600, height = 900)
print(pp_all)
dev.off()


