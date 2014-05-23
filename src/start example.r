# oposSOM start template
#
# set preferences
p <- opossom.new(list(dataset.name = "Unnamed",
                      error.model = "all.samples",
                      dim.1stLvlSom = 20,
                      dim.2ndLvlSom = 20,
                      training.extension = 1,
                      rotate.SOM.portraits = 0,
                      flip.SOM.portraits = F,
                      database.dataset = "",
                      database.id.type = "",
                      geneset.analysis = F,
                      geneset.analysis.exact = F,
                      max.parallel.cores = detectCores() / 2,
                      spot.threshold.samples = 0.65,
                      spot.coresize.modules = 3,
                      spot.threshold.modules = 0.95,
                      spot.coresize.groupmap = 5,
                      spot.threshold.groupmap = 0.75,
                      feature.centralization = T,
                      sample.quantile.normalization = T,
                      pairwise.comparison.list = list())

# load Data
p$indata <- .......

# define some groups
p$group.labels <- .....

# define some sample colors
# pipeline$group.colors
#   <- c("col1","col2",...)[match(p$group.labels, unique(p$group.labels))]

# execute
opossom.run(p)
