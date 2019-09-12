#Plot Hydrophobicity ----
library(reshape2)
library(ggplot2)

df <- data.frame('Mean Free'<- grade_1_hydrophob$Total, 'G1_Total'<- pangle.data$mean_of_complex_groups)
head(df)



x10plusrange=c(0,70)
y10plusrange=c(0,12)

plot(x10plusrange,
     y10plusrange,
     main = 'Range of P-angles for 4plus(green) mapped with 10plus(brown) redundant Ab groups', 
     xlab ='All antibody groups with redundant structures', 
     ylab = 'Packing Angle',
     breaks = 20, 
     col='white',
     freq = FALSE
)

lines(range_of_all_groups, col='purple2', lwd=2)
lines(range_spaced_10plus_redundgroups, col='green', lwd=3)