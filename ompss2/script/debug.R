library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(scales)
library(jsonlite)
library(readr)
library(viridis)

dpi = 300
h = 3
w = 4

# ---------------------------------------------------------------------

df = read_delim("energy.csv", delim=",") %>%
  mutate(Ekindiff = Ekin - lag(Ekin)) %>%
  mutate(Epotdiff = Epot - lag(Epot)) %>%
  mutate(iter = iter + 1)

p = ggplot(df, aes(x=iter)) +
  geom_line(aes(y=Ekin + Epot, color="Total")) +
  geom_line(aes(y=Ekin, color="Kinetic")) +
  geom_line(aes(y=Epot, color="Potential")) +
#  geom_line(aes(y=Epotg, color="Potential ghosts")) +
  theme_bw() +
  labs(x="Timestep", y="Energy", title="Energy conservation", linetype="Energy type")

ggsave("energy.png", plot=p, width=w, height=h, dpi=dpi)

dfref = read_delim("ref-energy.csv", delim=",")

p = ggplot(df, aes(x=iter)) +
  geom_line(            aes(y=Etot, linetype="New",       color="Total")) +
  geom_line(            aes(y=Ekin, linetype="New",       color="Kinetic")) +
  geom_line(            aes(y=Epot, linetype="New",       color="Potential")) +
  geom_line(data=dfref, aes(y=Etot, linetype="Reference", color="Total")) +
  geom_line(data=dfref, aes(y=Ekin, linetype="Reference", color="Kinetic")) +
  geom_line(data=dfref, aes(y=Epot, linetype="Reference", color="Potential")) +
#  geom_line(aes(y=Epotg, color="Potential ghosts")) +
  theme_bw() +
  labs(x="Timestep", y="Energy", title="Energy conservation",
    color="Energy Type", linetype="Variant")

ggsave("energy-ref.png", plot=p, width=w, height=h, dpi=dpi)

p = ggplot(df, aes(x=iter)) +
  geom_line(aes(y=df$Etot - dfref$Etot)) +
  theme_bw() +
  labs(x="Timestep", y="Energy", title="Energy error")

ggsave("energy-error.png", plot=p, width=w, height=h, dpi=dpi)

p = ggplot(df, aes(x=iter)) +
  geom_line(aes(y=Etot - Etot[1], linetype="Corrected")) +
  geom_line(data=dfref, aes(y=Etot - Etot[1], linetype="Reference")) +
  geom_hline(yintercept = 0.0005 * df$Etot[1], color = "red") +
  geom_hline(yintercept = -0.0005 * df$Etot[1], color = "red") +
  theme_bw() +
  labs(x="Timestep", y="Energy", title="Total energy")

ggsave("etot.png", plot=p, width=w, height=h, dpi=dpi)

p = ggplot(df, aes(x=iter)) +
  geom_line(aes(y=Epot)) +
  theme_bw() +
  labs(x="Timestep", y="Energy", title="Potential energy")

ggsave("epot.png", plot=p, width=w, height=h, dpi=dpi)

p = ggplot(df, aes(x=iter)) +
  geom_line(aes(y=Ekin)) +
  theme_bw() +
  labs(x="Timestep", y="Energy", title="Kinetic energy")

ggsave("ekin.png", plot=p, width=w, height=h, dpi=dpi)

p = ggplot(df, aes(x=iter)) +
  geom_line(aes(y=Ekindiff, color="Kinetic")) +
  geom_line(aes(y=Epotdiff, color="Potential")) +
  theme_bw() +
  labs(x="Timestep", y="Diff energy", title="Diff energy")

ggsave("ekin-delta.png", plot=p, width=w, height=h, dpi=dpi)

dfd = read_delim("energy-drift.csv", delim=",")

p = ggplot(df, aes(x=iter)) +
  geom_line(data=dfd, aes(y=Etot - Etot[1], color="Original")) +
  geom_line(data=df, aes(y=Etot - Etot[1], color="Corrected")) +
  theme_bw() +
  labs(x="Timestep", y="Delta energy", title="Total energy", color="Version")

ggsave("etot-drift.png", plot=p, width=w, height=h, dpi=dpi)

# ---------------------------------------------------------------------

df = read_delim("dhist.csv", delim=",") %>%
  mutate(lcount = log(count + 1)) %>%
  mutate(capcount = ifelse(count > 0, 1, 0))

p <- ggplot(df, aes(x = iter+1, y=limit)) +
  geom_tile(aes(fill = lcount)) +
  theme_bw() +
  scale_fill_viridis() +
  labs(fill = "Log count") +
  geom_hline(yintercept = 2.5, color = "red") +
  geom_hline(yintercept = 2.5+0.3, color = "blue") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(x="Iteration", y="Distance", title="Distance histogram in box 0")

ggsave("dhist.png", plot=p, width=w, height=h, dpi=dpi)

# ---------------------------------------------------------------------

df = read_delim("fhist.csv", delim=",") %>%
  mutate(lcount = log(count + 1))

p <- ggplot(df, aes(x = iter + 1, y=limit)) +
  geom_tile(aes(fill = lcount)) +
  theme_bw() +
  scale_fill_viridis() +
  labs(fill = "Log count") +
  labs(x="Iteration", y="Force", title="Force histogram in box 0")

ggsave("fhist.png", plot=p, width=w, height=h, dpi=dpi)

# ---------------------------------------------------------------------

df = read_delim("vhist.csv", delim=",") %>%
  mutate(lcount = log(count + 1))

p <- ggplot(df, aes(x = iter + 1, y=limit)) +
  geom_tile(aes(fill = lcount)) +
  theme_bw() +
  scale_fill_viridis() +
  labs(fill = "Log count") +
  labs(x="Iteration", y="Velocity", title="Velocity histogram in box 0")

ggsave("vhist.png", plot=p, width=w, height=h, dpi=dpi)

## ---------------------------------------------------------------------
#
N = 20
df = read_delim("atompos.csv", delim=",") %>%
  filter(atom < N) %>%
#  filter(ghost == 0) %>%
  mutate(ghost = as.factor(ghost)) %>%
  mutate(atom = as.factor(atom))


L = 6.7183847655300291
R = 2.8
k = 0.05
p = ggplot(df, aes(x=x, y=y, color=iter, shape=ghost)) +
  geom_hline(yintercept = 0, color = "blue") +
  geom_hline(yintercept = L, color = "blue") +
  geom_vline(xintercept = 0, color = "blue") +
  geom_vline(xintercept = L, color = "blue") +
  geom_hline(yintercept = 0+R, color = "gray") +
  geom_hline(yintercept = L-R, color = "gray") +
  geom_vline(xintercept = 0+R, color = "gray") +
  geom_vline(xintercept = L-R, color = "gray") +
  geom_point() +
#  geom_line(aes(group=atom)) +
  scale_color_viridis() +
  coord_fixed() +
  theme_bw() +
  labs(x="x", y="y", title=sprintf("Movement of %d atoms in XY", N))

ggsave("atompos-xy.png", plot=p, width=w, height=h, dpi=dpi)

p = ggplot(df, aes(x=z, y=y, color=iter, shape=ghost)) +
  geom_hline(yintercept = 0, color = "blue") +
  geom_hline(yintercept = L, color = "blue") +
  geom_vline(xintercept = 0, color = "blue") +
  geom_vline(xintercept = L, color = "blue") +
  geom_hline(yintercept = 0+R, color = "gray") +
  geom_hline(yintercept = L-R, color = "gray") +
  geom_vline(xintercept = 0+R, color = "gray") +
  geom_vline(xintercept = L-R, color = "gray") +
  geom_point() +
#  geom_line(aes(group=atom)) +
  scale_color_viridis() +
  coord_fixed() +
  theme_bw() +
  labs(x="z", y="y", title=sprintf("Movement of %d atoms in YZ", N))

p = ggplot(df, aes(x=x, y=z, color=iter, shape=ghost)) +
  geom_hline(yintercept = 0, color = "blue") +
  geom_hline(yintercept = L, color = "blue") +
  geom_vline(xintercept = 0, color = "blue") +
  geom_vline(xintercept = L, color = "blue") +
  geom_hline(yintercept = 0+R, color = "gray") +
  geom_hline(yintercept = L-R, color = "gray") +
  geom_vline(xintercept = 0+R, color = "gray") +
  geom_vline(xintercept = L-R, color = "gray") +
  geom_point() +
#  geom_line(aes(group=atom)) +
  scale_color_viridis() +
  coord_fixed() +
  theme_bw() +
  labs(x="x", y="z", title=sprintf("Movement of %d atoms in XZ", N))

ggsave("atompos-xz.png", plot=p, width=w, height=h, dpi=dpi)

#
#p = ggplot(df, aes(x=iter+1, y=neigh, color=atom)) +
#  geom_line() +
#  theme_bw() +
#  labs(x="Timestep", y="Number of nearby atoms", title="Nearby atoms")
#
#ggsave("neigh.png", plot=p, width=w, height=h, dpi=dpi)
