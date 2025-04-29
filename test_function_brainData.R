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

dir.out <- 'C:/project/power analysis/FG_EM/'

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
    M_candidates = 2:10,
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
  
  ## output results
  saveRDS(rst, file = paste0(dir.out, 'humanbrain_151673_', REGIONS[i], '_FG_EM_fit.rds'))
  
}

# Print summary
cat("\nTiming Summary:\n")
print(timing_results)

# Optionally save results
write.csv(timing_results, paste0(dir.out,"model_selection_timing_results_humanBrain151673.csv"), row.names = FALSE)
plot_density_FGM(300, rst3$best_model)





########running time for original FG 

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
  FGM <- FG_mixture(dat.sub[,1:2],Iter = 20000)
  
  # End timing
  time_end <- Sys.time()
  
  # Calculate duration
  duration <- as.numeric(time_end - time_start, units = "secs")
  
  # Store results
  timing_results$time_seconds[i] <- duration
  
  
  # Print progress
  cat(sprintf("Region: %s, Time: %.2f seconds",
              REGIONS[i], duration))
  
}

# Print summary
cat("\nTiming Summary:\n")
print(timing_results)

# Optionally save results
write.csv(timing_results, paste0(dir.out,"FG_timing_results_humanBrain151673.csv"), row.names = FALSE)









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


dir.out <- 'C:/project/power analysis/FG_EM/'

### WM
REGIONS <- sort(unique(dat$region))
regime <- REGIONS
regime[2] <- 'compactLV'
regime[5] <-'TrabecularLV'

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
    M_candidates = 2:7,
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
  
  ## output results
  saveRDS(rst, file = paste0(dir.out, 'chickenHeart_', regime[i], '_FG_EM_fit.rds'))
  
}

# Print summary
cat("\nTiming Summary:\n")
print(timing_results)

# Optionally save results
write.csv(timing_results, paste0(dir.out,"model_selection_timing_results_chickenHeart_M2_7.csv"), row.names = FALSE)



rst <- FG_EM(x = as.matrix(region.norm[, 1:2]), M = 3, max_iter = 1000, tol = 1e-3)


rst3 <- select_best_M(x = as.matrix(dat.sub[, 1:2]), M_candidates = 2:6, max_iter = 1000, tol = 1e-1)




########### chicken heart

dir.in <- 'C:/Users/xie15/Downloads/FG_EM_simulation/'

