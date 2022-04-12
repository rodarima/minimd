library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(scales)
library(jsonlite)
library(readr)
library(tidyr)
library(viridis)

dpi = 200
h = 4
w = 12

# ---------------------------------------------------------------------

#dfbb = read_delim("counters/bb.trace", delim=",") %>%
#  mutate(type="Blocking") %>% 
#  mutate(time = time * 1e-9) %>%
#  mutate(reltime = time - time[1])

#dfnn = read_delim("counters/nn-4nodes.trace", delim=",") %>%
dfnn = read_delim("counters/nn-4nodes-noreneigh.trace", delim=",") %>%
  mutate(type="Non-blocking") %>%
  mutate(time = time * 1e-9) %>%
  mutate(reltime = time - time[1]) %>%
  filter(reltime > 16) %>%
  filter(reltime < 17)
#  filter(reltime > 20.6) %>%
#  filter(reltime < 21.7)

df = bind_rows(dfnn) %>%
  mutate(as.factor(type)) %>%
  mutate(deltatime = reltime - lag(reltime)) %>%
  mutate(bwsend = deltasend / deltatime / (1024 * 1024)) %>%
  mutate(bwrecv = deltarecv / deltatime / (1024 * 1024)) %>%
  pivot_longer(c(bwsend, bwrecv), names_to="direction", values_to="bw")

print(df)

p = ggplot(df, aes(x=reltime, color=direction)) +
  geom_line(aes(y=bw)) +
  geom_point(aes(y=bw)) +
  facet_grid(direction ~ .) +
  theme_bw() +
  #scale_y_continuous(trans = 'log10')+
  labs(x="Time", y="MiB/s", title="Bandwidth", color="Traffic")

ggsave("net-4nodes.png", plot=p, width=w, height=h, dpi=dpi)
