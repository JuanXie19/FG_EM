library(Seurat)
library(igraph)
library(ggplot2)
library(dplyr)

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


DAT <- as.data.frame(t(counts))
DAT$region <- brain$benmarklabel

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
dat.sub <- dat %>% filter(region =='WM')


rst <- FG_EM(x = as.matrix(dat.sub[,1:2]), M = 4, max_iter =1000, tol = 1e-4 )

rst2 <- FG_EM(x = as.matrix(dat.sub[,1:2]), M = 5, max_iter =1000, tol = 1e-4 )


##### subset first, and then normalize
dat <- data.frame(row = coords[,1],col = coords[,2],region = brain$benmarklabel)
dat.sub <- dat %>% filter(region =='WM')

region.norm <- norm_coords(as.matrix(dat.sub[,1:2]),xmin = 0, xmax = 1,ymin = 0,ymax = 1)
rst <- FG_EM(x = as.matrix(region.norm[, 1:2]), M = 5, max_iter = 1000, tol = 1e-4)
rst <- FG_EM(x = as.matrix(region.norm[, 1:2]), M = 3, max_iter = 1000, tol = 1e-3)


rst3 <- select_best_M(x = as.matrix(region.norm[, 1:2]), M_candidates = 2:6, max_iter = 1000, tol = 1e-4)

plot_density_FGM(300, rst3$best_model)