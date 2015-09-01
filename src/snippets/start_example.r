
library(oposSOM)

# oposSOM start template
#
# Set preferences

env <- opossom.new(list(dataset.name = "Unnamed",
                      error.model = "all.samples.LPE",

                      dim.1stLvlSom = 20,
                      dim.2ndLvlSom = 20,

                      training.extension = 1,
                      rotate.SOM.portraits = 0,
                      flip.SOM.portraits = F,

                      database.dataset = "auto",
                      database.id.type = "auto",

                      geneset.analysis = T,
                      geneset.analysis.exact = F,
                      geneset.analysis.samplespots = F,

                      spot.threshold.samples = 0.65,
                      spot.coresize.modules = 3,
                      spot.threshold.modules = 0.95,
                      spot.coresize.groupmap = 5,
                      spot.threshold.groupmap = 0.75,

                      feature.centralization = T,
                      sample.quantile.normalization = T,

                      pairwise.comparison.list = list() ) )


# Load input data
env$indata <- .......

# Define sample groups
env$group.labels <- .....

# Define sample colors

	# env$group.colors
	#   <- c("col1","col2",...)[match(env$group.labels, unique(env$group.labels))]

# execute

opossom.run(env)










#### Tissue Example from Vignette ####

library(oposSOM)

env <- opossom.new(list(dataset.name="Tissues", dim.1stLvlSom=20))
data(opossom.tissues)

env$indata <- opossom.tissues

env$group.labels <- c(rep("Homeostasis", 2),"Endocrine","Digestion","Exocrine","Epithelium",
                          "Reproduction","Muscle",rep("Immune System", 2),rep("Nervous System", 2) )
env$group.colors <- c(rep("gold", 2),"red2","brown","purple","cyan","pink","green2",rep("blue2", 2),rep("gray", 2) )

opossom.run(env)



