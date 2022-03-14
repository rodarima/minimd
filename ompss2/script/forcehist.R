library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(scales)
library(jsonlite)
library(readr)
library(viridis)

# Load the arguments (argv)
args = commandArgs(trailingOnly=TRUE)
if (length(args)>0) { input_file = args[1] } else { input_file = "forcehist.csv" }
if (length(args)>1) { output = args[2] } else { output = "?" }

df = read_delim(input_file, delim=",", comment="#") %>%
  mutate(lcount = log(count + 1)) %>%
  mutate(pcount = log((count) * (bin + 1)))
#  filter(bin > 20)

print(df)

dpi = 300
h = 4
w = 6

# ---------------------------------------------------------------------

p <- ggplot(df, aes(x = iter, bin)) +
  geom_tile(aes(fill = lcount)) +
  theme_bw() +
  scale_fill_viridis() +
  labs(fill = "Log count") +
#  geom_vline(xintercept = 29.5, color = "red", size=0.1) +
  labs(x="Timestep", y="Force bin", title="Force histogram in box 0")

ggsave(sprintf("%s.heatmap.png", input_file), plot=p, width=w, height=h, dpi=dpi)
