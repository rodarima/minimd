library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(scales)
library(jsonlite)
library(readr)

# Load the arguments (argv)
args = commandArgs(trailingOnly=TRUE)
if (length(args)>0) { input_file = args[1] } else { input_file = "energy.csv" }
if (length(args)>1) { output = args[2] } else { output = "?" }

df = read_delim(input_file, delim=",")

print(df)

dpi = 300
h = 4
w = 6

# ---------------------------------------------------------------------

p = ggplot(df, aes(x=iter+1)) +
  geom_line(aes(y=Ekin + Epot*3, linetype="Total"), color="red") +
  geom_line(aes(y=Ekin, linetype="Kinetic")) +
  geom_line(aes(y=Epot*3, linetype="Potential")) +
  theme_bw() +
  labs(x="Timestep", y="Energy", title="Energy conservation", linetype="Energy type")

ggsave(sprintf("%s.png", input_file), plot=p, width=w, height=h, dpi=dpi)
