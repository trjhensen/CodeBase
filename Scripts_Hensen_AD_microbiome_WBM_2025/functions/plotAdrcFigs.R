# Plot ADRC figures


# - Define a helper function 'using' to load and install packages if needed
using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}

setwd('C://Users//mspg//Documents//ADRC')

# - Load required packages
using('tidyverse','here','pheatmap','RColorBrewer','viridis','tidytext','ggalluvial','ggplotify')


################################################################################
# Display heatmap of flux-microbe correlations

# Load dataset
data <- read.csv(
  here('outputs','microbetoflux','enetFluxMicrobeRes.csv'),
  check.names = FALSE,row.names = 1) %>%
  rownames_to_column(var = "Species") %>% 
  mutate(Species = str_replace(Species,'_',' ')) %>%
  replace(is.na(.), 0)

data <- data %>% filter(Species != 'Flux-associated taxa')
data <- data %>% filter(Species != 'Sex')

# Remove formaldehyde
# data <- data %>% select(-Formaldehyde)

if (F) { # Do not apply now
  # Filter on microbes that correlate above the threshold
  threshold <- 0.1
  fun <- function(x,t) {abs(data[sapply(x, is.numeric)]) > t }
  data <- data[apply(fun(data,threshold), 1, any),]
}

# Split correlation data and annotations
dataProcessed <- as.data.frame(data[,colnames(data)])
rownames(dataProcessed) <- data$Species
dataProcessed <- dataProcessed %>% select(-Species)


# Create pheatmap
#spec <- annotationProcessed %>% select(-Species) #%>% relocate(,'Phylum',.after = 'Genus')
#rownames(spec) <- annotationProcessed$Species
#spec <- spec %>% select(-Order,-Family,-Genus)# %>% relocate(,'Phylum',.after = 'Family')

numCols <- 9
# mycol <- c(brewer.pal(numCols,"Blues")[numCols:1],"white",brewer.pal(numCols,"Reds")[1:numCols])
# breaks <- rev(c(0.8, 0.7,0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0, -0.1, -0.2, -0.3, -0.4,-0.5,-0.6,-0.7,-0.8,-0.9))+0.05
mycol <- c(brewer.pal(numCols,"Blues")[3:2],"white","white",brewer.pal(numCols,"Oranges")[2:numCols-1])
breaks <- rev(c(0.8,0.7,0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0, -0.1, -0.2,-0.3,-0.4))+0.1#+0.05

#dataProcessed[abs(dataProcessed)<0.009] <- 0
dataProcessed <- round(dataProcessed, digits = 2)

# Create heatmap and save figure
#png(filePath, width = 5, height = 6, units = 'in', res = 450)
#png(filePath, width = 6, height = 8, units = 'in', res = 450)
plt <- pheatmap(dataProcessed,
                #color = colorRampPalette(c("blue", "white", "red"))(10),
                color = mycol,
                breaks = breaks,
                display_numbers = T,
                fontsize_number = 7,
                number_color = "black",
                number_format = "%0.2g",
                cluster_rows = T, cluster_cols = F,
                clustering_method="average",
                clustering_distance_rows='euclidean',
                clustering_distance_cols='euclidean',         
                fontsize_row=9,
                fontsize_col=9,
                fontsize = 8,
                #annotation_row = spec,
                #annotation_colors = ann_colors,
                angle_col = 315,
                #main = "Correlations of species relative abundances\n with metabolite blood fluxes in mmol/day/person",
)

ggplt <- as.ggplot(plt) +
  labs(title = "Standardised elastic net coefficients for the\npredicted microbe-metabolite associations") +
  #guides(fill = guide_legend(ncol=1)) +
  theme(plot.title = element_text(size = rel(1)))

print(ggplt)
ggsave(here('outputs','figures','fluxMicrobeENET.png'), device='png',  width=5, height = 6, units = 'in')


#numCols <- 10
#mycol <- c(brewer.pal(numCols,"Blues")[2:1],"white",brewer.pal(numCols,"Reds")[1:(numCols-2)])
#breaks <- rev(c(0.9,0.8,0.7,0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0, -0.1, -0.2,-0.3))+0.05
#breaks <- rev(c(0.8,0.6, 0.4, 0.2, 0,  -0.2))+0.05

#breaks <- rev(c(0, -0.1, -0.2, -0.3, -0.4,-0.5,-0.6,-0.7,-0.8,-0.9))+0.05
# ann_colors = list(
#   `Shift in PD` = c(`Higher in PD`="#440154", `Lower in PD`="#21908C",Unchanged="#FDE725"),
#   Phylum = c(Actinobacteria="#fcfdbf", Bacteroidetes="#fe9f6d", Euryarchaeota="#de4968",
#              Firmicutes="#8c2981",Proteobacteria="#3b0f70",Verrucomicrobia="#000004"))


# Create heatmap and save figure
filePath <- here('outputs','figures','fluxMicrobeENET.png')
#png(filePath, width = 5, height = 6, units = 'in', res = 450)
png(filePath, width = 10.5, height = 4, units = 'in', res = 400)
plt <- pheatmap(t(dataProcessed),
                color = mycol,
                breaks = breaks,
                display_numbers = T,
                number_color = 'black',
                fontsize_number = 10,
                cluster_rows = F, cluster_cols = T,
                clustering_method="average",
                clustering_distance_rows='euclidean',
                clustering_distance_cols='euclidean',         
                fontsize_row=10,
                fontsize_col=10,
                fontsize = 11,
                angle_col = "315",
                #annotation_row = spec,
                #annotation_colors = ann_colors,
                main = "Elastic net standardised coefficients for microbial relative abundances\nagainst the predicted fluxes"
)
print(plt)
dev.off()     


# ################################################################################
# # Produce bar plots for the mean flux sensitivity upon microbial abundances
# 
# # Load dataset
# data <- read.csv(
#   here('outputs','microbetoflux','microbeSensitivityTable.csv'),
#   check.names = FALSE) %>%
#   # Process data and make the negative values as positive (prediction reduction in flux)
#   mutate(Taxa = str_replace(Taxa,'_',' ')) %>%
#   mutate(
#     Metabolite = factor(Metabolite,levels = c('S-Adenosyl-L-methionine','L-arginine','Creatine','Taurine','Formate')))
# 
# 
# # Filter on the top n microbes per metabolite
# n <- 10
# dataPruned <- data %>% group_by(Metabolite) %>% 
#   arrange(desc(Mean),.by_group = TRUE) %>%
#   slice(1:n) %>% 
#   mutate(Taxa = factor(Taxa,levels = Taxa))
# 
# # Create bar plots
# 
# # Start r graphics png
# fileName <- here('outputs','figures','fluxSensitivityPlot.png')
# png(file=fileName,width = 15, height = 5, units = "in", res = 400)
# 
# # Plot
# ggplot(dataPruned) +
#   geom_col( aes(x=Mean, y=reorder_within(Taxa,Mean,Metabolite)),fill='darkorange2',colour='black') +
#   geom_errorbarh(aes(y = reorder_within(Taxa,Mean,Metabolite), xmax = `97.5CI`, xmin = `2.5CI`,),height = .3) +
#   scale_y_reordered() +
#   labs(y = '', 
#        x = '',#Predicted flux contributions in mmol/day/person',
#        title = 'Microbial contributions towards predicted fluxes in blood',
#        subtitle = 'Predicted maximal flux contributions in mmol/day/person',#paste0('Top microbial contributors'),
#        fill = 'Change in relative abundance\nin PD microbiomes') +
#   facet_wrap(vars(Metabolite),scales = "free",ncol=4,labeller = "label_both") + 
#   scale_fill_viridis(discrete = TRUE, option='viridis') +
#   theme_bw()+
#   theme(legend.position="bottom",
#         axis.text.y = element_text(colour='black'))
# 
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# ################################################################################
# # Create chord diagram of top microbe-metabolite interactions
# 
# topClusterCombinations <- read.csv(
#   here('outputs','microbetoflux','mostImportantMicrobialLimiters.csv'),
#   check.names=F) %>% 
#   rename(Species = Row) %>%
#   pivot_longer(-Species, names_to = 'Metabolite') %>%
#   mutate(Species = gsub('_',' ',Species),
#          Metabolite = gsub('_',' ',Metabolite)) %>%
#   filter(value!=0) 
# 
# # Start r graphics png
# fileName <- here('outputs','figures','topMicrobeMets_alluvium.png')
# png(file=fileName,width = 10, height = 6, units = "in", res = 450)
# 
# nameSize <- 4.6
# ggplot(data = topClusterCombinations,
#        aes(axis1 = Species, axis2 = Metabolite, y=value)) +
#   geom_alluvium(aes(fill = Metabolite)) +
#   geom_stratum(fill='grey90') +
#   geom_text(stat = "stratum", 
#             aes(label = after_stat(stratum)),size=nameSize) +
#   scale_x_discrete(limits = c("Species", "Metabolite"),expand = c(0, 0),
#                    labels=c("Species" = "Gut microbial species", 
#                             "Metabolite" = "Predicted metabolite\nin blood")) +
#   scale_fill_viridis_d() +
#   labs(title = 'Predicted key influencers of metabolic productions in blood',y='') +
#   theme_minimal() +
#   theme(#legend.position = "none",
#     plot.title = element_text(size=16),
#     axis.text.y=element_blank(),
#     axis.text.x=element_text(size=14),
#     panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank(),
#     legend.position = 'bottom',
#     legend.text = element_text(size=14),
#     legend.title = element_text(size=14))
# dev.off()
