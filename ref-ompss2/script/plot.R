library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(scales)
library(jsonlite)
library(readr)

# Load the arguments (argv)
args = commandArgs(trailingOnly=TRUE)
if (length(args)>0) { input_file = args[1] } else { input_file = "out.csv" }
if (length(args)>1) { output = args[2] } else { output = "?" }

df = read_delim(input_file, delim=" ", comment="#") %>%
  mutate(K = T * 3/2)

print(df)

dpi = 300
h = 4
w = 6

# ---------------------------------------------------------------------

p = ggplot(df, aes(x=Step)) +
  geom_line(aes(y=K,    linetype="Kinetic energy (K)")) +
  geom_line(aes(y=U,    linetype="Potential energy (U)")) +
  geom_line(aes(y=K+U,  linetype="Total energy (K + U)")) +
  theme_bw() +
  labs(x="Timestep", y="Energy per particle U/N", title="Energy conservation in miniMD")

ggsave(sprintf("%s.energy.png", input_file), plot=p, width=w, height=h, dpi=dpi)

# ---------------------------------------------------------------------

p = ggplot(df, aes(x=Step)) +
  geom_line(aes(y=K+U, color="Total")) +
  theme_bw() +
  labs(x="Timestep", y="Energy per particle U/N", title="Energy conservation in miniMD")

ggsave(sprintf("%s.totenergy.png", input_file), plot=p, width=w, height=h, dpi=dpi)

# ---------------------------------------------------------------------

p = ggplot(df, aes(x=Step)) +
  geom_line(aes(y=Vcenter, color="Vc")) +
  theme_bw()

ggsave(sprintf("%s.vc.png", input_file), plot=p, width=w, height=h, dpi=dpi)
