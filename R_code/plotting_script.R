#Archive --------------------
library(reshape2)
library(ggplot2)

df <- data.frame('Mean Free'<- pangle.data$mean_of_free_groups, 'Mean Complex'<- pangle.data$mean_of_complex_groups)
head(df)                  

compare_packing_angles <- 
  melt(pangle.data$mean_of_free_groups, 
       pangle.data$mean_of_complex_groups)

ggplot(compare_packing_angles, 
       aes(compare_packing_angles, fill = Packing_Angles)) +
      geom_density(alpha = 0.2)



#Dataframe ------------------- ALL

pangle.data <- 
  data.frame(all_ab_groups, mean_of_all_groups, mean_of_free_groups, mean_of_complex_groups,
             range_of_all_groups, range_of_free_groups, range_of_complex_groups,
             sd_of_all_groups, sd_of_free_groups, sd_of_complex_groups)


all_ab_groups=seq(1,66)

head(pangle.data)


  

#10plus redundancies range-------------------

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


#10plus redundancies SD-------------------

x10plussd=c(0,70)
y10plussd=c(0,6)

plot(x10plussd,
     y10plussd,
     main = 'SD of pangles for 4plus(green) mapped with 10plus(brown) redundant Ab groups', 
     xlab ='All antibody groups with redundant structures', 
     ylab = 'Packing Angle',
     breaks = 20, 
     col='white',
     freq = FALSE
)

lines(sd_of_all_groups, col='gold', lwd=2)
lines(sd_spaced_10plus_redundgroups, col='deeppink', lwd=3)



#10plus redundancies SD vs Range----------


x10plussdrange=c(0,70)
y10plussdrange=c(0,12)

plot(x10plussdrange,
     y10plussdrange,
     main = 'SD and Range across 
     4+ and 10+ redundancy Ab groups', 
     xlab ='Antibody groups', 
     ylab = 'Packing Angle',
     breaks = 20, 
     col='white',
     freq = FALSE
)
legend(10, 12, legend=c("Range_4plus", "Range_10plus",'SD_4plus','SD_10plus'),
       col=c("purple2", "green",'gold','deeppink'), lty=1:2, cex=0.8)


lines(sd_of_all_groups, col='gold', lwd=2)
lines(sd_spaced_10plus_redundgroups, col='deeppink', lwd=3)

lines(range_of_all_groups, col='purple2', lwd=2)
lines(range_spaced_10plus_redundgroups, col='green', lwd=3)

#Range only ---------------------

  #All groups range
ggplot(data=pangle.data, aes(x=All_Ab_groups,
                             y=range_of_all_groups,
                             group=66)) +
  geom_line(col='purple')+
  geom_line(y=range_of_complex_groups, col='red')+
  geom_line(y=range_of_free_groups, col='blue')+
  geom_point(col='purple')+
  ggtitle('Range of P-Angles by redundant Antibody groups')+
  xlab("All Antibody groups")+
  ylab('Range in Degrees')

  #Comp vs Free range linegraph

ggplot(data=pangle.data, aes(x=All_Ab_groups,
                             y=range_of_complex_groups,
                             group=66)) +
  geom_line(col='red')+
  geom_line(y=range_of_free_groups, col='blue')+
  ggtitle('Range of P-Angles by redundant Antibody groups')+
  xlab("All Antibody groups")+
  ylab('Range in Degrees')

  #comp vs Free range histogram

hist(pangle.data$range_of_all_groups, 
     main='Range of P-angle for Ab groups', 
     xlab ='Packing Angle', col='purple', 
     xlim=c(0, 10), breaks = 20, freq = FALSE)+
    hist(pangle.data$range_of_free_groups, 
     xlab ='Packing Angle', breaks = 20, 
     freq = FALSE)+
    hist(pangle.data$range_of_complex_groups, 
     xlab ='Packing Angle', breaks = 20, 
     freq = FALSE)

melt_all<-melt(pangle.data$mean_of_free_groups)

compare_pangle_range <- 
  melt(pangle.data$range_of_free_groups,
        pangle.data$range_of_complex_groups)

melted_groups<-seq(1,594)

ggplot(melt_all, aes(x=All_Ab_groups,
                     y=pangle.data$range_of_complex_groups,
                     fill=All_Ab_groups
                     ))+
  geom_line(col=('red'))+
  geom_area(col='red', alpha=0.5, fill='red')+
  geom_line(y=pangle.data$range_of_free_groups, col='blue')+
  ggtitle('Range of P-Angles by Free and Complexed Antibodies')+
    xlab("All Antibody groups")+
    ylab('Range in Degrees')
  
legend(1, 10, legend=c("Complex",'Free'),
         col=c("red", "blue"), lty=1:2, cex=0.8)

legend(10, 12, legend=c("Range_4plus", "Range_10plus"),
       col=c("purple2", "green"), lty=1:2, cex=0.8)
  
  
  geom_area(y=pangle.data$mean_of_complex_groups, col='red')

         
show(all_ab_groups)
#-----------------------
xrange=c(0,70)
yrange=c(0,12)

plot(xrange,
     yrange,
     main = 'Range of all Packing angles by Ab groups', 
     xlab ='All redundant antibody structure groups', 
     ylab = 'Range of Packing Angles',
     breaks = 20, 
     col='white',
     freq = FALSE
)

lines(range_of_complex, col='red')
lines(range_of_free, col='blue')
lines(range_of_all, col='purple')


boxplot(range_of_all, 
        main = 'Range of All pdb Packing angles', 
        xlab ='Packing Angles', 
        notch=TRUE,
        col='purple3',
        freq = FALSE)


hist(range_of_complex, 
     main = 'Range of all Complex pdb Packing angles', 
     xlab ='Packing Angle', 
     xlim=c(0,10), 
     breaks = 50, 
     col='purple3',
     freq = FALSE)



#Mean only ---------------------

  #All groups mean scatter DONE
ggplot(data=pangle.data, aes(x=All_Ab_groups,
                             y=pangle.data$mean_of_all_groups,
                             group=66)) +
  geom_line(col='purple')+
  geom_point(col='purple')+
  ggtitle('Mean of P-Angles by redundant Antibody groups')+
  xlab("All Antibody groups")+
  ylab('Mean in Degrees')


  #All groups mean histogram DONE
theme_set(theme_bw())
ggplot(pangle.data, aes(x=mean_of_all_groups, y=All_Ab_groups), 
       col='purple', 
       breaks=100) +
  geom_count(col='purple')+
  ggtitle('Mean of P-angles by redundant Antibody groups')+
  xlab('Mean of All Ab groups')+
  ylab('Frequency')


  #Complex vs Free Mean density------------
mean_all<- pangle.data$mean_of_all_groups
mean_free <- pangle.data$mean_of_free_groups
mean_complex <-pangle.data$mean_of_complex_groups

#mean_free$pangle<- 'Free'
#mean_complex$pangle <- 'Complex'

pangledegrees <- rbind(mean_free, mean_complex)

ggplot(melt_all, aes(mean_all, fill=pangledegrees))+
                    geom_density(alpha=0.5, col='red', fill='red')
      geom_density(aes(mean_complex, alpha=0.3))

  #Complex vs Free Mean scatter------------

ggplot(melt_all, aes(x=All_Ab_groups,
                     y=pangle.data$mean_of_complex_groups
))+
  geom_line(col=('red'))+
  geom_line(y=pangle.data$mean_of_free_groups, col='blue')+
  ggtitle('Mean of free vs complex Packing angles by Ab groups')+
  xlab("All Antibody groups")+
  ylab('Mean of P-angles in Degrees')

quantile(mean_of_free_groups)

hist(mean_of_free_groups,
     main = 'Mean across 67 Free Antibody groups', 
     xlab ='Packing Angle', 
     breaks = 60, 
     col='orangered',
     freq = FALSE
)

#SD variation----------------

xsd<-c(0,66)
ysd<-c(0,5)

plot(xrange,
     yrange,
     main = 'SD/Variation across 66 Free vs Complex Antibody groups', 
     xlab ='Free (blue) vs Complex(red) Antibody groups', 
     ylab = 'SD/Variation of Packing Angle',
     breaks = 20, 
     col='white',
     freq = FALSE
)

lines(sd_of_free_groups, col='blue')
lines(sd_of_complex_groups, col='red')
lines(sd_of_all_groups, col='purple')

#Archive-----------------------

