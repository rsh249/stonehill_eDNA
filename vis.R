library(ggplot2)
library(stringr)
library(taxonomizr)
library(dplyr)
library(ggsci)

##Script for processing Kraken2 read classifications
# see kraken_script for commands to run kraken2



getNamesAndNodes()
# getAccession2taxid(types=c('nucl_gb'))
# getAccession2taxid()
# system("gunzip *.gz")
# read.accession2taxid(list.files('.','accession2taxid'),'accessionTaxa.sql')
# print(paste('taxonomizr database built and located at', getwd(), sep=' '))


taxaNodes<-read.nodes.sql('./nodes.dmp') #takes time?
taxaNames<-read.names.sql('./names.dmp')

files = list.files('~/ITS2_meta/', pattern='.tab', full.names = T)

rbuild = data.frame(barcode=character(), percent=numeric(), 
                    reads=numeric(), 
                    readsUnique=numeric(), 
                    rank=factor(), taxonID=numeric(), 
                    superkingdom=character(), 
                    phylum=character(), class=character(), 
                    order=character(), 
                    family=character(), genus=character(), 
                    species=character(),
                    perctotal=numeric())
rsum = data.frame(genus=character(), readsUnique=numeric())
for (f in 1:length(files)){
  print(f)
  r=read.csv(files[f], sep='\t', stringsAsFactors = F, header=F)
  rID = getId(str_trim(r$V6), sqlFile = 'nameNode.sqlite')
  
  rGen = getTaxonomy(
    rID,
    sqlFile = 'nameNode.sqlite',
    desiredTaxa = c(
      "superkingdom",
      "phylum",
      "class",
      "order",
      "family",
      "genus",
      "species"
    )
  )
  r = cbind(rep(paste('barcode_0', f, sep=''), nrow(r)), r,rID, rGen)
  colnames(r) = c('barcode', 'percent', 
                  'reads', 
                  'readsUnique', 
                  'rank', 'taxonID', 'name', 'ID2', 'superkingdom', 
                  'phylum', 'class', 'order', 
                  'family', 'genus', 'species')
  sum.plants = sum(filter(r, phylum=='Streptophyta') %>% select(readsUnique))
  rf = filter(r, phylum=='Streptophyta') %>% 
    filter(genus != "NA") %>%
    group_by(genus) %>%
    summarize(readsSum = sum(readsUnique)) %>%
    filter(readsSum >= 0.05*sum.plants)
  
  rf = cbind(rep(paste('barcode_0', f, sep =''), nrow(rf)), rf)
  colnames(rf) = c('barcode', 'genus', 'readsUnique')
  

  rbuild = rbind(rbuild, r)
  rsum = rbind(rsum, rf)
  
}


options(scipen=999, digits=1)

ggplot(data = rsum, aes(x=barcode, 
                        y=readsUnique, 
                        fill = genus)
       ) +
  geom_bar(position="fill", stat="identity") +
  scale_color_npg() + 
  scale_fill_npg() +
  theme_linedraw() + 
  theme(
    axis.text = element_text(size = 8),
    axis.text.x  = element_text(angle = 70, hjust=1),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  ) + xlab('') + ylab('Relative Proportion of Plant Reads')

ggsave('Figure2.png', height=4, width = 7.25)


ggplot(data = rbuild) + 
  geom_col(aes(x=barcode, y=readsUnique)) +
  theme_linedraw() +
  ylab('Plant Classified Reads') + 
  xlab('')

ggsave('Figure1.png', height=4, width = 7.25)


