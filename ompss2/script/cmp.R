library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(scales)
library(jsonlite)
library(readr)

# Load the arguments (argv)
args = commandArgs(trailingOnly=TRUE)

input_file1 = args[1]
input_file2 = args[2]

df1 = read_delim(input_file1, delim=" ", comment="#") # %>% mutate(T=1.5*T)
df2 = read_delim(input_file2, delim=" ", comment="#") # %>% mutate(T=1.5*T)

dpi = 300
h = 4
w = 6

# ---------------------------------------------------------------------

p = ggplot() +
  geom_line(data=df1, aes(x=Step, y=T, color="Thermal (T)", linetype=input_file1)) +
  geom_line(data=df1, aes(x=Step, y=U, color="Potential (U)", linetype=input_file1)) +
  geom_line(data=df1, aes(x=Step, y=T+U, color="Total (T+U)", linetype=input_file1)) +
  geom_line(data=df2, aes(x=Step, y=T, color="Thermal (T)", linetype=input_file2)) +
  geom_line(data=df2, aes(x=Step, y=U, color="Potential (U)", linetype=input_file2)) +
  geom_line(data=df2, aes(x=Step, y=T+U, color="Total (T+U)", linetype=input_file2)) +
  labs(x="Timestep", y="Energy", title="Energy conservation in miniMD")

ggsave(sprintf("%s.cmp.energy.png", input_file1), plot=p, width=w, height=h, dpi=dpi)

