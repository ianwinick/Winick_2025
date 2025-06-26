library(tidyverse)

################################################################################
# This code calculates and plots the beta deviations for all severity-year     #
# combinations using the bootstrapped beta diversity values and the randomized #
# (null) beta diversity values. The beta deviations are calculated by          #
# subtracting the mean null beta diversity value from each bootstrapped beta   #
# diversity value and dividing by the standard deviation of the null beta      #
# diversity.                                                                   #
#                                                                              #  
# Regression analysis is at the end of the code, including model selection     #
# (using AIC) to select a linear versus quadratic model.                       #
################################################################################

# Read in bootstrapped observed beta diversity values and null beta diversity values
obs <- read_csv("Bootstrap Observed Beta Values/bootstrap_beta.csv")

null <- read_csv("Null Model/nullBeta_allyears.csv") %>%
  pivot_longer(cols=c(tax, fun), names_to="type", values_to="null") %>%
  group_by(year, severity, type) %>%
  summarise(sd=sd(null), null=mean(null)) %>%
  select(year, type, severity, null, sd)

# Join observed and null beta diversity values and calculate beta deviation
data <- obs %>%
  left_join(null) %>%
  mutate(dev=NA, across(dev, ~(beta-null)/sd)) %>%
  mutate(severity=factor(severity, levels=c("U", "L", "H")))

# Calculate mean and standard deviation of beta deviations for use in plots
deviation <- data %>%
  group_by(year, type, severity) %>%
  summarise(sd=sd(dev), dev=mean(dev))

# Plot beta deviations through time by severity
tiff("Fig1.tiff", width=18, height=18, units="cm", res=600)

deviation %>%
  mutate(severity=case_when(
    severity=="U" ~ "Unburned",
    severity=="L" ~ "Low",
    severity=="H" ~ "High"
  )) %>%
  mutate(severity=factor(severity, levels=c("Unburned", "Low", "High"))) %>%
  mutate(type=case_when(
    type=="tax" ~ "Taxonomic",
    type=="fun" ~ "Functional"
  )) %>%
  ggplot(aes(year, dev, color=severity)) +
  geom_hline(yintercept=0, color="gray50", linetype="dashed") +
  geom_point(size=2, show.legend=F) +
  geom_errorbar(aes(ymin=dev-sd, ymax=dev+sd), width=0.05, show.legend=F) +
  geom_smooth(aes(linetype=severity), method="lm", formula=y~poly(x,2), se=FALSE) +
  facet_wrap(~factor(type, levels=c("Taxonomic", "Functional")), ncol=2) +
  scale_linetype_manual(values=c("solid", "dotdash", "dotted")) + 
  scale_color_manual(values=c("cadetblue3", "darkgoldenrod3", "black")) +
  scale_y_continuous(breaks=seq(-40,10,2),
                     labels = c(-40, rep("",4),
                                -30, rep("",4),
                                -20, rep("",4),
                                -10, rep("",4),
                                0, rep("",4), 10)) +
  xlab("Year since fire") +
  ylab("β-deviation") +
  labs(color="Fire severity", linetype="Fire severity") +
  theme_classic() +
  theme(axis.text.x=element_text(color="black", size=8, family="Arial"), 
        axis.text.y=element_text(color="black", size=8, family="Arial"),
        axis.title.x=element_text(color="black", size=10, family="Arial"),
        axis.title.y=element_text(color="black", size=10, family="Arial"),
        strip.text.x = element_text(color="black", size = 10, family="Arial"),
        legend.text=element_text(color="black", size=8, family="Arial"),
        legend.title=element_text(color="black", size=10, family="Arial"),
        legend.key.width=unit(1.2,"cm"),
        legend.position="bottom")

dev.off()

# Linear and quadratic regressions, with AIC test
# Based on AIC test, quadratic regression is used in final results
lin_tax <- lm(dev~year*severity, data[data$type=="tax",])
lin_fun <- lm(dev~year*severity, data[data$type=="fun",])
qua_tax <- lm(dev~year*severity*I(year^2), data[data$type=="tax",])
qua_fun <- lm(dev~year*severity*I(year^2), data[data$type=="fun",])

summary(lin_tax)
summary(lin_fun)
summary(qua_tax)
summary(qua_fun)

AIC(lin_tax)
AIC(lin_fun)
AIC(qua_tax)
AIC(qua_fun)
