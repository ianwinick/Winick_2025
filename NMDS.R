library(tidyverse)
library(vegan)
library(FD)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# Read in 2023 community and trait data
data <- read_csv("CommunityMatrix2023.csv") %>%
  column_to_rownames(var="plot")

traits <- read_csv("traitMatrix.csv") %>%
  filter(spp %in% colnames(data)) %>%
  column_to_rownames(var="spp")

# Create vector of fire severities for use in PERMANOVA groups
severity <- data %>%
  mutate(severity = case_when(
    row.names(.) %in% 1:20 ~ "Unburned",
    row.names(.) %in% 21:40 ~ "Low",
    row.names(.) %in% 41:60 ~ "High"
  )) %>%
  select(severity)

################################################################################
# Taxonomic composition ########################################################
################################################################################
tax_nmds <- metaMDS(data, distance="bray", k=2, trymax=100)
tax_nmds$stress
stressplot(tax_nmds)

tax_scores <- as_tibble(scores(tax_nmds, display="sites"), rownames="sites") %>%
  mutate(severity=case_when(
    sites %in% c(1:20) ~ "Unburned",
    sites %in% c(21:40) ~ "Low",
    sites %in% c(41:60) ~ "High"
  ))

# Calculate group centroids
tax_avg <- tax_scores %>%
  mutate(severity=factor(severity, levels=c("Unburned", "Low", "High"))) %>%
  group_by(severity) %>%
  dplyr::summarize(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))

tax_scores %>%
  mutate(severity=factor(severity, levels=c("Unburned", "Low", "High"))) %>%
  ggplot(aes(x=NMDS1, y=NMDS2, linetype=severity, color=severity)) +
  stat_ellipse(geom="polygon", aes(group=severity, fill=severity), alpha = 0.3, size=.75) +
  geom_point(shape=1, size=2) +
  geom_point(data=tax_avg, aes(NMDS1, NMDS2, shape=severity), size=3, color="black") +
  scale_linetype_manual(values=c("solid", "dotdash", "dotted")) +
  scale_color_manual(values=c("cadetblue3", "darkgoldenrod3", "black")) +
  scale_fill_manual(values=c("cadetblue3", "darkgoldenrod3", "black")) +
  theme_classic()

# Taxonomic PERMANOVA
tax_permanova <- adonis2(vegdist(data, method="bray") ~ severity, permutations=9999, data=severity)
pairwise.adonis(vegdist(data, method="bray"), severity$severity, perm=9999)

# Taxonomic test of beta dispersion
anova(betadisper(vegdist(data, method="bray"), severity$severity, type="centroid"))
TukeyHSD(betadisper(vegdist(data, method="bray"), severity$severity, type="centroid"))

################################################################################
# Community weighted means for use in functional composition ###################
################################################################################

# Log transform traits
traitMatrix <- traits %>%
  mutate(ldmc=log(ldmc)) %>%
  mutate(height=log(height)) %>%
  mutate(sla=log(sla))

# Visualize distribution of transformed trait values
hist(traitMatrix$ldmc)
hist(traitMatrix$height)
hist(traitMatrix$sla)

# Calculate community weighted means of traits
cwm <- dbFD(traitMatrix, data)$CWM

################################################################################
# Functional composition #######################################################
################################################################################
trait_nmds <- metaMDS(cwm, distance="euclidean", k=2, trymax=100)
trait_nmds$stress
stressplot(trait_nmds)

fun_scores <- as_tibble(scores(trait_nmds, display="sites"), rownames="sites") %>%
  mutate(severity=case_when(
    sites %in% c(1:20) ~ "Unburned",
    sites %in% c(21:40) ~ "Low",
    sites %in% c(41:60) ~ "High"
  ))

# Calculate group centroids
fun_avg <- fun_scores %>%
  group_by(severity) %>%
  summarise(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))

fun_scores %>%
  ggplot(aes(x=NMDS1, y=NMDS2, linetype=severity, color=severity)) +
  stat_ellipse(geom="polygon", aes(group=severity, fill=severity), alpha = 0.3, size=.75) +
  geom_point(shape=1, size=2) +
  geom_point(data=fun_avg, aes(NMDS1, NMDS2, shape=severity), size=3, color="black") +
  scale_linetype_manual(values=c("solid", "dotdash", "dotted")) +
  scale_color_manual(values=c("cadetblue3", "darkgoldenrod3", "black")) +
  scale_fill_manual(values=c("cadetblue3", "darkgoldenrod3", "black")) +
  labs(color="", fill="", linetype="", shape="") +
  theme_classic()

# Functional PERMANOVA
fun_permanova <- adonis2(vegdist(cwm, method="euclidean") ~ severity, permutations=9999, data=severity)
pairwise.adonis(vegdist(cwm, method="euclidean"), severity$severity, perm=9999)

# Functional test of beta dispersion
anova(betadisper(vegdist(cwm, method="euclidean"), severity$severity, type="centroid"))
TukeyHSD(betadisper(vegdist(cwm, method="euclidean"), severity$severity, type="centroid"))

################################################################################
# Bind taxonomic and functional NMDS scores and print faceted NMDS plot ########
################################################################################
tax_scores <- tax_scores %>% mutate(type="tax")
fun_scores <- fun_scores %>% mutate(type="fun")
scores <- rbind(tax_scores, fun_scores)
tax_avg <- tax_avg %>% mutate(type="tax")
fun_avg <- fun_avg %>% mutate(type="fun")
avg <- rbind(tax_avg, fun_avg)

tiff("Fig3.tiff", width=8.5, height=14, units="cm", res=600)

scores %>%
  mutate(severity=factor(severity, levels=c("Unburned", "Low", "High"))) %>%
  ggplot(aes(x=NMDS1, y=NMDS2, linetype=severity, color=severity)) +
  stat_ellipse(geom="polygon", aes(group=severity, fill=severity), alpha = 0.3, size=.75) +
  geom_point(shape=1, size=2) +
  geom_point(data=avg, aes(NMDS1, NMDS2, shape=severity), size=3, color="black") +
  scale_linetype_manual(values=c("solid", "dotdash", "dotted")) +
  scale_color_manual(values=c("cadetblue3", "darkgoldenrod3", "black")) +
  scale_fill_manual(values=c("cadetblue3", "darkgoldenrod3", "black")) +
  facet_wrap(~factor(type, levels=c("tax", "fun")), ncol=1, scales="free") +
  labs(color="", fill="", linetype="", shape="") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=8, family="Arial"), 
        axis.text.y=element_text(color="black", size=8, family="Arial"),
        axis.title.x=element_text(color="black", size=10, family="Arial"),
        axis.title.y=element_text(color="black", size=10, family="Arial"),
        legend.text=element_text(color="black", size=8, family="Arial"),
        legend.title=element_text(color="black", size=10, family="Arial"),
        legend.position="bottom")

dev.off()