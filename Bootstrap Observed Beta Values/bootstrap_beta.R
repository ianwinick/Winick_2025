library(tidyverse)

boot2020 <- read_csv("Bootstrap Observed Beta Values/bootstrap_2020.csv")
boot2021 <- read_csv("Bootstrap Observed Beta Values/bootstrap_2021.csv")
boot2022 <- read_csv("Bootstrap Observed Beta Values/bootstrap_2022.csv")
boot2023 <- read_csv("Bootstrap Observed Beta Values/bootstrap_2023.csv")

# Bind all bootstrap data frames together (2020-2023)
data <- rbind(boot2020, boot2021, boot2022, boot2023) %>%
  mutate(severity=factor(severity, levels=c("U", "L", "H")))

# Export bound data frame
write.csv(data, "Bootstrap Observed Beta Values/bootstrap_beta.csv", row.names=FALSE)

################################################################################
# Data visualization ###########################################################
################################################################################

# Calculate mean bootstrapped beta diversity and standard deviation
data <- data %>%
  group_by(year, type, severity) %>%
  summarize(sd=sd(beta), beta=mean(beta))

# Plot it. Note, these are the observed beta diversity values and are not yet compared to the null (i.e. random) expectation.
data[data$type=="tax",] %>%
  ggplot(aes(year, beta)) +
  geom_point() +
  geom_errorbar(aes(ymin=beta-sd, ymax=beta+sd), width=0.1) +
  facet_wrap(~severity) +
  ylab("mean Bray-Curtis dissimilarity") +
  xlab("year since fire") +
  theme_bw()

data[data$type=="fun",] %>%
  ggplot(aes(year, beta)) +
  geom_point() +
  geom_errorbar(aes(ymin=beta-sd, ymax=beta+sd), width=0.1) +
  facet_wrap(~severity) +
  ylab("functional beta diversity (Q)") +
  xlab("year since fire") +
  theme_bw()
