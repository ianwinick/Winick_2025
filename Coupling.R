library(tidyverse)
library(vegan)
library(FD)
library(extrafont)

# R code for measuring the coupling between taxonomic and function diversity. Coupling is measured using functional redundancy, represented as the slope between species richness and Rao's quadratic entropy (Q).

# Read in community matrices and trait matrix
df20 <- read_csv("CommunityMatrix2020.csv") %>%
  filter(!plot %in% c(33,42,45,46,48,51,53,56,57)) %>%
  mutate(year=1)

df21 <- read_csv("CommunityMatrix2021.csv") %>%
  mutate(year=2)

df22 <- read_csv("CommunityMatrix2022.csv") %>%
  mutate(year=3)

df23 <- read_csv("CommunityMatrix2023.csv") %>%
  mutate(year=4)

df24 <- read_csv("CommunityMatrix2024.csv") %>%
  mutate(year=5)

traitMatrix <- read_csv("traitMatrix.csv") %>%
  column_to_rownames(var="spp") %>%
  arrange(row.names(.))

# Bind community matrices across years, attach severity classes
library(plyr)
df <- tibble(rbind.fill(df20, df21, df22, df23, df24)) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  select(sort(colnames(.))) %>%
  mutate(severity=factor(case_when(
    plot %in% 1:20 ~ "U",
    plot %in% 21:40 ~ "L",
    plot %in% 41:60 ~ "H"
  ), levels=c("U", "L", "H"))) %>%
  select(plot, year, severity, everything(.))
detach("package:plyr", unload = TRUE)

# Calculate Simpson's diversity
simp <- diversity(df[,4:69], index="simpson")

# Calculate Q
Q <- dbFD(traitMatrix, df[,4:69], scale.RaoQ=T, message=F)$RaoQ

# Combine richness and Q in one data frame
fRed <- tibble("plot"=df$plot, "year"=df$year, "simp"=simp, "Q"=Q) %>%
  mutate(severity=factor(case_when(
    plot %in% 1:20 ~ "U",
    plot %in% 21:40 ~ "L",
    plot %in% 41:60 ~ "H"
  ), levels=c("U", "L", "H"))) %>%
  unite(sev_year, c(severity, year), sep="_") %>%
  filter(Q>0)

# Calculate regression values for each severity-year
sev_year <- c("U_1", "U_2", "U_3", "U_4", "U_5", "L_1", "L_2", "L_3", "L_4", "L_5", "H_1", "H_2", "H_3", "H_4", "H_5")
reg.e = NULL
reg.r2 = NULL
reg.ar2 = NULL
reg.p = NULL

for(i in sev_year){
  x <- fRed %>%
    #mutate(Q=(Q-mean(Q))/sd(Q), rich=(rich-mean(rich))/sd(rich)) %>%
    filter(sev_year==i)
  
  reg <- summary(lm(Q~simp, x))
  reg.e[i] = coef(reg)["simp", "Estimate"]
  reg.r2[i] = reg$r.squared
  reg.ar2[i] = reg$adj.r.squared
  reg.p[i] = coef(reg)["simp", "Pr(>|t|)"]
}

# Extract regression values and put in a data frame
reg.data <- data.frame("e"=reg.e, "r2"=reg.r2, "ar2"=reg.ar2, "p"=reg.p)
reg.data <- rownames_to_column(reg.data, var="id")
reg.data <- reg.data %>%
  separate(id, c("severity", "year")) %>%
  mutate(severity=factor(severity, levels=c("U", "L", "H")))
write.csv(reg.data, "coupling_data.csv")

reg.data %>% filter(p<0.051)

reg.data %>% slice_min(e)
reg.data %>% slice_max(e)

reg.data %>%
  group_by(severity) %>%
  reframe(cv=sd(e)/mean(e))

reg.data %>%
  select(severity, year, e) %>%
  pivot_wider(names_from=severity, values_from=e)

# Functional redundancy (Q~species richness) plotted by severity
tiff("Fig2.tiff", width=18, height=14, units="cm", res=600)

fRed %>%
  separate(sev_year, into=c("severity", "year")) %>%
  mutate(severity=case_when(
    severity=="U" ~ "Unburned",
    severity=="L" ~ "Low",
    severity=="H" ~ "High"
  ))%>%
  mutate(severity=factor(severity, levels=c("Unburned", "Low", "High"))) %>%
  mutate(year=factor(year)) %>%
  ggplot(aes(simp, Q, color=year, linetype=year)) +
  geom_point(show.legend=F, color="gray80") +
  geom_smooth(method="lm", se=F) +
  facet_wrap(~severity, ncol=3) +
  scale_color_manual(values=c("darkolivegreen4", "firebrick3", "deepskyblue3", "darkgoldenrod3", "black")) +
  scale_linetype_manual(values=c("solid", "longdash", "dashed", "dotdash", "dotted")) +
  xlab("Simpson's Diversity Index") +
  ylab("Rao's quadratic entropy (Q)") +
  labs(color="Year since fire", linetype="Year since fire") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=8, family="Arial"), 
        axis.text.y=element_text(color="black", size=8, family="Arial"),
        axis.title.x=element_text(color="black", size=10, family="Arial"),
        axis.title.y=element_text(color="black", size=10, family="Arial"),
        strip.text.x = element_text(size = 10, family="Arial"),
        legend.key.width=unit(1.2,"cm"),
        legend.position="bottom")
  
 dev.off()
 