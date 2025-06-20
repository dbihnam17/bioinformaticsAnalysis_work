### This script was testing out creating fishplots from our WGS data (without pyclone)

## I am using the fishplot library (https://github.com/chrisamiller/fishplot)
## with the code shown on the github page to model our data

# Load fishplot library
library(fishplot)

# Identify timepoints of samples     
timepoints=c(1, 2)      

# Plug in the clonal VAFs per clone and timepoint
frac.table = matrix(
  c(12, 46, 0,
    0, 49, 12),
  ncol=length(timepoints))

# Assign parents to each clone
parents = c(0, 0, 2)

# CreatE a Fish object
fish = createFishObject(frac.table,parents,timepoints=timepoints,
                        clone.labels=c("Founding", "Subclone 1","Subclone 2"), 
                        fix.missing.clones=TRUE)

fish = layoutClones(fish)


# The sample times you want to be plotted
sample.times = c(1, 2)
# The sample time labels
vlabel = c("Baseline", "Progression")

# Function to plot the Fish Plot
fishPlot(fish,shape="spline",title.btm="sample_ID",
         vlines=sample.times, vlab=vlabel, cex.title=0.5, bg.type="solid", col.vline = "#2D6D66")
