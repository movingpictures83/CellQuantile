# combine k-mer correlation tests, filter for significant values and species resolution
library(dplyr)
library(ggplot2)
library(scales)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")
source("RIO.R")

input <- function(inputfile) {
   pfix <- prefix()
   parameters <<- readParameters(inputfile)
   celllines <<- parameters["celllines", 2]
   samplename <<- parameters["samplename", 2]
   myc <<- readRDS(paste(pfix, parameters["c", 2], sep="/"))
   c2 <<- readRDS(paste(pfix, parameters["c2", 2], sep="/"))
   kr <<- readRDS(paste(pfix, parameters["kr", 2], sep="/"))
   qtile <<- as.numeric(parameters["qtile", 2])
}

run <- function() {
c3 <<- left_join(myc, c2, by = 'name') %>% 
  left_join(select(kr, rank, name) %>% distinct()) %>% 
  subset(r1>0 & r2>0 & r3>0 & p<0.05 & p1<0.05 & p2<0.05 & p3<0.05 & rank == 'S')
cell.lines = read.csv(celllines, header=T, sep=' ')
df = cell.lines[,1:10] %>% mutate(study = 'cell lines'); #df = df[, -2]
df = rbind(df, kr)

ggplot(subset(df, name %in% c3$name), aes(rpmm, fill = study, ..scaled..)) + 
         geom_density(alpha = 0.5, color = NA) +
         facet_grid(~name, scales='free') + 
         scale_fill_manual(values = c(`cell lines` = 'grey50', zhang = 'navyblue')) +
         theme_minimal() + 
         xlab('Microbiome reads per million') +
         ylab('Density') +
         scale_x_log10(labels = trans_format("log10", math_format(10^.x)),
                       breaks = trans_breaks("log10", function(x) 10^x, n=4),
                       oob = scales::squish, expand = c(0,0)) +
         scale_y_continuous(expand = c(0,0)) +
         theme(legend.title = element_blank(),
               panel.border = element_rect(fill = NA, color = 'black'),
               strip.background = element_blank(),
               axis.ticks.x = element_line(size=0),
               axis.ticks.y = element_blank(),
               panel.grid = element_blank(), 
               strip.text = element_text(color = 'black', size = 10),
               axis.text.x = element_text(color = 'black', size = 10),
               axis.text.y = element_text(color = 'black', size = 10),
               axis.title.y = element_text(color = 'black', size = 10),
               axis.title.x = element_text(color = 'black', size = 10),
               plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
               legend.key.size = unit(0.2, "cm"),
               legend.text = element_text(color = 'black', size = 10),
               legend.position = 'bottom')

q_df = cell.lines %>%
  group_by(name, rank) %>%
  summarize(CLrpmm = 10^quantile(log10(rpmm), qtile, na.rm = T),
            .groups = 'keep')

left_join(c3, q_df, by = c('name', 'rank')) %>%
  left_join(subset(kr, sample == samplename) %>% select(name, rpmm), by = c('name'))
}


output <- function(outputfile) {
write.table(c3$taxid, paste(outputfile, "taxa.tsv", sep="/"), row.names=F, col.names=F)
saveRDS(c3, paste(outputfile, "c3.rds", sep="/"))
}
