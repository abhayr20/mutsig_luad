#Load libraries
library(tidyverse)
library(dplyr)

setwd("/data/abhay/young_onset/")

stats = read.delim('signatures/SPA/Assignment_Solution/Solution_Stats/Assignment_Solution_Samples_Stats.txt')

ggplot(stats) +
  aes(x=Cosine.Similarity) +
  labs(x='')+
  geom_histogram(aes(y = after_stat(density))) +
  geom_density(col = 4, lwd = 1.5) +
  geom_vline(aes(xintercept = 0.9),
             col = 2, lwd = 1.5) +
  labs(x = 'Cosine Similarity') +
  theme_bw()

# Read activities matrix
acts = read.delim('signatures/SPE/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt')
avg_acts = colMeans(acts[,-1])
barplot(avg_acts)

# Generate stacked barplot (percent stacked)
acts_tidy = acts %>%
  pivot_longer(cols = !Samples,
               names_to = 'Signature',
               values_to = 'Mutations')

ggplot(acts_tidy) +
  aes(x = Samples, y = Mutations, fill = Signature) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

number_of_mutations = rowSums(acts[,-1])

# Selecting the activities of only the top 10 mutated cases
top_10_mutated_samples = acts[order(number_of_mutations,
                                    decreasing = T)[1:10],]

# Reformatting and plotting
top_10_mutated_samples %>%
  pivot_longer(cols = !Samples,
               names_to = 'Signature',
               values_to = 'Mutations') %>%
  ggplot() +
  aes(x = reorder(Samples, Mutations), y = Mutations, fill = Signature) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  labs(x = 'Samples')  +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
