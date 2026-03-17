library(dplyr)
library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)
library(gridExtra)


BAT <- read.csv("BAT_MEs_matchedmice.csv", header = TRUE, stringsAsFactors = TRUE)
Liver <- read.csv("Liver_MEs_matchedmice.csv", header = TRUE, stringsAsFactors = TRUE)
Muscle <- read.csv("Muscle_MEs_matchedmice.csv", header = TRUE, stringsAsFactors = TRUE)

liver1 = ggplot(Liver, aes(x=factor(Diet, level=c("CTL", "ValR")), MEgreen, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "Liver.MEgreen", x = "Diet", y = "")

bat1 = ggplot(BAT, aes(x=factor(Diet, level=c("CTL", "ValR")), MEred, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "BAT.MEred", x = "Diet", y = "")

bat2 = ggplot(BAT, aes(x=factor(Diet, level=c("CTL", "ValR")), MEblue, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "BAT.MEblue", x = "Diet", y = "")

liver2 = ggplot(Liver, aes(x=factor(Diet, level=c("CTL", "ValR")), MEtan, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "Liver.MEtan", x = "Diet", y = "")

muscle1 = ggplot(Muscle, aes(x=factor(Diet, level=c("CTL", "ValR")), MEturquoise, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "Muscle.MEturquoise", x = "Diet", y = "")

muscle2 = ggplot(Muscle, aes(x=factor(Diet, level=c("CTL", "ValR")), MEred, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "Muscle.MEred", x = "Diet", y = "")

liver3 = ggplot(Liver, aes(x=factor(Diet, level=c("CTL", "ValR")), MEbrown, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "Liver.MEbrown", x = "Diet", y = "")

liver4 = ggplot(Liver, aes(x=factor(Diet, level=c("CTL", "ValR")), MEred, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "Liver.MEredn", x = "Diet", y = "")

muscle3 = ggplot(Muscle, aes(x=factor(Diet, level=c("CTL", "ValR")), MEgreen, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "Muscle.MEgreen", x = "Diet", y = "")


plot_grid(liver1,bat1,bat2,liver2,muscle1,muscle2,liver3,liver4,muscle3,ncol=3)





















bat1 = ggplot(BAT, aes(x=factor(Diet, level=c("CTL", "ValR")), MEred, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "MEred", x = "Diet", y = "")

bat4 = ggplot(BAT, aes(x=factor(Diet, level=c("CTL", "ValR")), MEturquoise, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "MEturquoise", x = "Diet", y = "")

bat5 = ggplot(BAT, aes(x=factor(Diet, level=c("CTL", "ValR")), MEsalmon, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "MEsalmon", x = "Diet", y = "")

bat6 = ggplot(BAT, aes(x=factor(Diet, level=c("CTL", "ValR")), MEmidnightblue, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "MEmidnightblue", x = "Diet", y = "")

plot_grid(bat1, bat2, bat3, bat4, bat5, bat6, ncol = 2)



liver1 = ggplot(Liver, aes(x=factor(Diet, level=c("CTL", "ValR")), MEgreen, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "MEgreen", x = "Diet", y = "")

liver2 = ggplot(Liver, aes(x=factor(Diet, level=c("CTL", "ValR")), MEpink, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "MEpink", x = "Diet", y = "")

liver3 = ggplot(Liver, aes(x=factor(Diet, level=c("CTL", "ValR")), MEtan, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "MEtan", x = "Diet", y = "")

liver4 = ggplot(Liver, aes(x=factor(Diet, level=c("CTL", "ValR")), MEblack, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "MEblack", x = "Diet", y = "")

liver5 = ggplot(Liver, aes(x=factor(Diet, level=c("CTL", "ValR")), MEroyalblue, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "MEroyalblue", x = "Diet", y = "")

liver6 = ggplot(Liver, aes(x=factor(Diet, level=c("CTL", "ValR")), MEcyan, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "MEcyan", x = "Diet", y = "")

plot_grid(liver1, liver2, liver3, liver4, liver5, liver6, ncol = 2)


muscle1 = ggplot(Muscle, aes(x=factor(Diet, level=c("CTL", "ValR")), MEtan, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "Muscle.MEtan", x = "Diet", y = "")

muscle2 = ggplot(Muscle, aes(x=factor(Diet, level=c("CTL", "ValR")), MEyellow, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "Muscle.MEyellow", x = "Diet", y = "")


plot_grid(muscle1,muscle2,ncol=2)



muscle3 = ggplot(Muscle, aes(x=factor(Diet, level=c("CTL", "ValR")), MEbrown, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "MEbrown", x = "Diet", y = "")

muscle4 = ggplot(Muscle, aes(x=factor(Diet, level=c("CTL", "ValR")), MEgrey60, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "MEgrey60", x = "Diet", y = "")

muscle5 = ggplot(Muscle, aes(x=factor(Diet, level=c("CTL", "ValR")), MEmagenta, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "MEmagenta", x = "Diet", y = "")

muscle6 = ggplot(Muscle, aes(x=factor(Diet, level=c("CTL", "ValR")), MEblack, fill = Diet)) + 
  geom_boxplot(position = position_dodge(width = 1)) + 
  geom_jitter(position = position_dodge(width = 1), alpha = 1) +
  facet_wrap(~Sex, nrow = 1) +
  theme(legend.position = "none") +
  labs(title = "MEblack", x = "Diet", y = "")

plot_grid(muscle1, muscle2, muscle3, muscle4, muscle5, muscle6, ncol = 2)




