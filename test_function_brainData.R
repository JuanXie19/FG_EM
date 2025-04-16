library(Seurat)
library(igraph)
library(ggplot2)
library(dplyr)
library(MASS)

dir.in <- 'C:/project/power analysis/151673/'

load(paste0(dir.in,'brain_processed.RData'))

meta <- read.csv(paste0(dir.in,'meta.csv'))

identical(colnames(brain),meta$barcode)

brain$benmarklabel <- meta$benmarklabel

regions <- sort(unique(brain$benmarklabel))
regions <- regions[nzchar(regions)]   # there are some unannotated spots, remove them

IND <- which(brain$benmarklabel %in% regions)
brain <- brain[,IND]

counts <- as.matrix(brain@assays$Spatial@counts)


## remove MT genes
mt_idx      <- grep("MT-",rownames(counts))
if(length(mt_idx)!=0){
  counts    <- counts[-mt_idx,]
}




coords <- brain@images$image@coordinates
coords <- coords[,1:2]

### normaliz first, and then subset

coords.norm <- norm_coords(as.matrix(coords[,1:2]),xmin = 0, xmax = 1,ymin = 0,ymax = 1)

dat <- data.frame(row = coords.norm[,1],col = coords.norm[,2],region = brain$benmarklabel)


###############
### step1: parameter estimation

## location

source('C:/Users/xie15/Downloads/CODE.R')

# scale the coordinates

ggplot(dat,aes(row,col,col=region))+geom_point()


### WM
REGIONS <- sort(unique(dat$region))

# Create a data frame to store timing results
timing_results <- data.frame(
  region = REGIONS,
  time_seconds = numeric(length(REGIONS)),
  best_M = numeric(length(REGIONS)),
  stringsAsFactors = FALSE
)

for (i in seq_along(REGIONS)) {
  dat.sub <- dat %>% filter(region == REGIONS[i])  # Fixed typo: 'regions' to 'region'
  
  # Start timing
  time_start <- Sys.time()
  
  # Run model selection
  rst <- select_best_M(
    x = as.matrix(dat.sub[, 1:2]),
    M_candidates = 2:6,
    max_iter = 1000,
    tol = 1e-1
  )
  
  # End timing
  time_end <- Sys.time()
  
  # Calculate duration
  duration <- as.numeric(time_end - time_start, units = "secs")
  
  # Store results
  timing_results$time_seconds[i] <- duration
  timing_results$best_M[i] <- rst$best_M_BIC  # or rst$best_M_AIC if you prefer
  
  # Print progress
  cat(sprintf("Region: %s, Time: %.2f seconds, Best M: %d\n",
              REGIONS[i], duration, timing_results$best_M[i]))
}

# Print summary
cat("\nTiming Summary:\n")
print(timing_results)

# Optionally save results
write.csv(timing_results, "model_selection_timing_results_humanBrain151673.csv", row.names = FALSE)





rst <- FG_EM(x = as.matrix(dat.sub[,1:2]), M = 4, max_iter =1000, tol = 1e-4 )




##### subset first, and then normalize
dat <- data.frame(row = coords[,1],col = coords[,2],region = brain$benmarklabel)
dat.sub <- dat %>% filter(region =='WM')
dat.sub <- dat %>% filter(region =='Layer4')



region.norm <- norm_coords(as.matrix(dat.sub[,1:2]),xmin = 0, xmax = 1,ymin = 0,ymax = 1)
rst <- FG_EM(x = as.matrix(region.norm[, 1:2]), M = 6, max_iter = 1000, tol = 1e-2)
rst <- FG_EM(x = as.matrix(region.norm[, 1:2]), M = 3, max_iter = 1000, tol = 1e-3)


rst3 <- select_best_M(x = as.matrix(region.norm[, 1:2]), M_candidates = 2:6, max_iter = 1000, tol = 1e-4)

plot_density_FGM(300, rst3$best_model)




############### test chicken heart data

dir.in <- 'C:/project/power analysis/chicken/'

load(paste0(dir.in,'heart_D14_seurat_processed.RData'))


regions <- sort(unique(heart_D14$region))

counts <- as.matrix(heart_D14@assays$Spatial@counts)


## remove MT genes
mt_idx      <- grep("MT-",rownames(counts))
if(length(mt_idx)!=0){
  counts    <- counts[-mt_idx,]
}



coords <- heart_D14@images$image@coordinates
coords <- coords[,1:2]

coords.norm <- norm_coords(as.matrix(coords[,1:2]),xmin = 0, xmax = 1,ymin = 0,ymax = 1)

dat <- data.frame(row = coords.norm[,1],col = coords.norm[,2],region = heart_D14$region)

reg.simp <- regions
reg.simp[2] <- 'Compact LV'
reg.simp[5] <- 'Trabecular Lv'


dat.sub <- dat %>% filter(region =='Atria') # work, M=3 is best when expon.scale =T, 
dat.sub <- dat %>% filter(region =='Compact LV and \ninter-ventricular septum') # work, M= 6
dat.sub <- dat %>% filter(region =='Epicardium') # work, M=2 is best, 
dat.sub <- dat %>% filter(region =='Right ventricle') # M=2 is best
dat.sub <- dat %>% filter(region =='Trabecular LV and \nendocardium') # 
dat.sub <- dat %>% filter(region =='Valves') # work, M=4

rst <- FG_EM(x = as.matrix(region.norm[, 1:2]), M = 3, max_iter = 1000, tol = 1e-3)


rst3 <- select_best_M(x = as.matrix(dat.sub[, 1:2]), M_candidates = 2:6, max_iter = 1000, tol = 1e-1)



#####  test running time
