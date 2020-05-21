
### run simulations
if (!file.exists("2d.rda"))
    source("2d.R")
if (!file.exists("lognormal_2d.rda"))
    source("lognormal_2d.R")
if (!file.exists("friedman.rda"))
    source("friedman.R")
if (!file.exists("lognormal_friedman.rda"))
    source("lognormal_friedman.R")
if (!file.exists("timings.rda"))
    source("timings.R")

### generate plots
pdf("plots.pdf", width = 6, height = 10)
source("plots.R")
dev.off()


