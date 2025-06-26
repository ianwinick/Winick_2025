library(tidyverse)

# This code pivots the community matrices into long-form, so they are formatted properly for the null model

# 2020 data
data2020 <- read_csv("CommunityMatrix2020.csv") %>%
 column_to_rownames(var="plot")

Community_longForm_2020 <- data2020 %>%
  rownames_to_column(var="plot") %>%
  pivot_longer(cols=-plot, names_to="spp", values_to="cover") %>%
  select(plot, spp, cover) %>%
  filter(cover>0)

# 2021 data
data2021 <- read_csv("CommunityMatrix2021.csv") %>%
  column_to_rownames(var="plot")

Community_longForm_2021 <- data2021 %>%
  rownames_to_column(var="plot") %>%
  pivot_longer(cols=-plot, names_to="spp", values_to="cover") %>%
  select(plot, spp, cover) %>%
  filter(cover>0)

# 2022 data
data2022 <- read_csv("CommunityMatrix2022.csv") %>%
  column_to_rownames(var="plot")

Community_longForm_2022 <- data2022 %>%
  rownames_to_column(var="plot") %>%
  pivot_longer(cols=-plot, names_to="spp", values_to="cover") %>%
  select(plot, spp, cover) %>%
  filter(cover>0)

# 2023 data
data2023 <- read_csv("CommunityMatrix2023.csv") %>%
  column_to_rownames(var="plot")

Community_longForm_2023 <- data2023 %>%
  rownames_to_column(var="plot") %>%
  pivot_longer(cols=-plot, names_to="spp", values_to="cover") %>%
  select(plot, spp, cover) %>%
  filter(cover>0)

write.csv(Community_longForm_2020, "Community_longForm_2020.csv", row.names=FALSE)
write.csv(Community_longForm_2021, "Community_longForm_2021.csv", row.names=FALSE)
write.csv(Community_longForm_2022, "Community_longForm_2022.csv", row.names=FALSE)
write.csv(Community_longForm_2023, "Community_longForm_2023.csv", row.names=FALSE)
