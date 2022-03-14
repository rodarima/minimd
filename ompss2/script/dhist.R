library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(scales)
library(jsonlite)
library(readr)
library(viridis)

df = read_delim("dhist.csv", delim=",", comment="#") %>%
  mutate(lcount = log(count + 1)) %>%
  mutate(pcount = log((count) * (bin + 1)))

dpi = 300
h = 4
w = 6

# ---------------------------------------------------------------------

p <- ggplot(df, aes(x = iter, y=limit)) +
  geom_tile(aes(fill = lcount)) +
  theme_bw() +
  scale_fill_viridis() +
  labs(fill = "Log count") +
  labs(x="Timestep", y="Distance", title="Distance histogram in box 0")

ggsave("dhist.png", plot=p, width=w, height=h, dpi=dpi)
