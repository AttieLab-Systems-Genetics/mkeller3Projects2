library(dplyr)
library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)
library(gridExtra)


data <- read.csv("Module_eigengenes_md20_annotated.csv", header = TRUE, stringsAsFactors = TRUE)


m1 = ggplot(data, aes(x=factor(Diet, level=c("Control", "Low AA", "Low ILE")), MEpurple, fill = Age)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  labs(title = "MEpurple", x = "Diet", y = "")

m2 = ggplot(data, aes(x=factor(Diet, level=c("Control", "Low AA", "Low ILE")), MEmagenta, fill = Age)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  labs(title = "MEmagenta", x = "Diet", y = "")

m3 = ggplot(data, aes(x=factor(Diet, level=c("Control", "Low AA", "Low ILE")), MEdarkgreen, fill = Age)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  labs(title = "MEdarkgreen", x = "Diet", y = "")

m4 = ggplot(data, aes(x=factor(Diet, level=c("Control", "Low AA", "Low ILE")), MEblack, fill = Age)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  labs(title = "MEblack", x = "Diet", y = "")

m5 = ggplot(data, aes(x=factor(Diet, level=c("Control", "Low AA", "Low ILE")), MEyellowgreen, fill = Age)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  labs(title = "MEyellowgreen", x = "Diet", y = "")

m6 = ggplot(data, aes(x=factor(Diet, level=c("Control", "Low AA", "Low ILE")), MEmidnightblue, fill = Age)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  labs(title = "MEmidnightblue", x = "Diet", y = "")

plot_grid(m1, m2, m3, m4, m5, m6, ncol = 3)

