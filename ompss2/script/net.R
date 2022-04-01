library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(scales)
library(jsonlite)
library(readr)
library(viridis)

dpi = 200
h = 4
w = 8

# ---------------------------------------------------------------------

dfbb = read_delim("counters/bb.trace", delim=",") %>%
  mutate(type="Blocking") %>% 
  mutate(time = time * 1e-9) %>%
  mutate(reltime = time - time[1])

dfnn = read_delim("counters/nn.trace", delim=",") %>%
  mutate(type="Non-blocking") %>%
  mutate(time = time * 1e-9) %>%
  mutate(reltime = time - time[1])

dfbbf = read_delim("counters/fixed-bb.trace", delim=",") %>%
  mutate(type="Blocking") %>% 
  mutate(time = time * 1e-9) %>%
  mutate(reltime = time - time[1]) %>%
  filter(reltime > 24.5) %>%
  mutate(reltime = time - time[1])

dfnnf = read_delim("counters/fixed-nn.trace", delim=",") %>%
  mutate(type="Non-blocking") %>%
  mutate(time = time * 1e-9) %>%
  mutate(reltime = time - time[1]) %>%
  filter(reltime > 6) %>%
  mutate(reltime = time - time[1])

dfpp = read_delim("counters/pingpong.trace", delim=",") %>%
  mutate(type="PingPong") %>% 
  mutate(time = time * 1e-9) %>%
  mutate(reltime = time - time[1])

df = bind_rows(dfbb, dfnn, dfpp) %>%
  mutate(as.factor(type)) %>%
  mutate(deltatime = reltime - lag(reltime)) %>%
  mutate(bwsend = deltasend / deltatime / (1024 * 1024)) %>%
  mutate(bwrecv = deltarecv / deltatime / (1024 * 1024))

print(df)

p = ggplot(df, aes(x=reltime)) +
  geom_line(aes(y=bwsend, color="Send")) +
  geom_line(aes(y=bwrecv, color="Recv")) +
  facet_grid(type ~ .) +
  theme_bw() +
  #scale_y_continuous(trans = 'log10')+
  labs(x="Time", y="MiB/s", title="Bandwidth", color="Traffic")

ggsave("net.png", plot=p, width=w, height=h, dpi=dpi)

df = bind_rows(dfbbf, dfnnf, dfpp) %>%
  mutate(as.factor(type)) %>%
  mutate(deltatime = reltime - lag(reltime)) %>%
  mutate(bwsend = deltasend / deltatime / (1024 * 1024)) %>%
  mutate(bwrecv = deltarecv / deltatime / (1024 * 1024))

print(df)

p = ggplot(df, aes(x=reltime, color=type)) +
  geom_line(aes(y=bwsend, color="Send")) +
  geom_line(aes(y=bwrecv, color="Recv")) +
  facet_grid(type ~ .) +
  theme_bw() +
  #scale_y_continuous(trans = 'log10')+
  labs(x="Time", y="MiB/s", title="Bandwidth (fixed)", color="Traffic")

ggsave("net-fixed.png", plot=p, width=w, height=h, dpi=dpi)

###############################################################

df2 = bind_rows(dfbbf, dfnnf) %>%
  mutate(as.factor(type)) %>%
  mutate(deltatime = reltime - lag(reltime)) %>%
  mutate(bwsend = deltasend / deltatime / (1024 * 1024)) %>%
  mutate(bwrecv = deltarecv / deltatime / (1024 * 1024))

p = ggplot(df2, aes(x=reltime, color=type)) +
  geom_line(aes(y=bwsend, color="Send")) +
  geom_line(aes(y=bwrecv, color="Recv")) +
  facet_grid(type ~ .) +
  theme_bw() +
  #scale_y_continuous(trans = 'log10')+
  labs(x="Time (s)", y="Bandwidth (MiB/s)", title="Bandwidth (fixed)", color="Traffic")

ggsave("net2.png", plot=p, width=w, height=h, dpi=dpi)
