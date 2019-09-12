#Tests:
install.packages("moments")
library(moments)
install.packages('nortest')
library(nortest)
library(dplyr)

#Mean Normalcy
kurtosis(mean_of_all_groups)
kurtosis(mean_of_free_groups)
kurtosis(mean_of_complex_groups)

skewness(mean_of_all_groups)
skewness(mean_of_free_groups)
skewness(mean_of_complex_groups)

#Redundancy tests
range_all<-pangle.data$range_of_all_groups

cor.test(All_redund_4_6, All_redund_7, method=c("pearson", "kendall", "spearman"))

wilcox.test(All_redund_4_6, All_redund_7, paired = TRUE)

ks.test(All_redund_4_6, All_redund_7)

var.test(All_redund_4_6, All_redund_7)

#Are the Ranges normally distributed
shapiro.test(All_redundRange_4_6)
shapiro.test(All_redundRange_7)
shapiro.test(range_of_free_groups)
shapiro.test(range_of_complex_groups)

t.test(h_chain, l_chain)
t.test(l_chain)

ks.test(range_of_complex_groups, range_of_free_groups)

skewness(All_redund_4_6)
kurtosis(All_redund_4_6)

skewness(All_redund_7)
kurtosis(All_redund_7)

skewness(range_of_free_groups)
kurtosis(range_of_free_groups)

skewness(range_of_complex_groups)
kurtosis(range_of_complex_groups)

ad.test(All_redund_4_6)
ad.test(All_redund_7)
ad.test(range_of_complex_groups)
ad.test(range_of_free_groups)

mean(range_of_all_groups)
mean(range_of_complex_groups)
mean(range_of_free_groups)

sd(range_of_complex_groups)
sd(range_of_free_groups)
sd(range_of_all_groups)

t.test(All_redund_4_6, All_redund_7)
t.test(range_of_complex_groups, range_of_free_groups)
#Are the ranges correlative/ variance test

var.test(All_redund_4_6, All_redund_7)
var.test(range_of_free_groups, range_of_complex_groups)

var.test(mean_of_complex_groups, mean_of_free_groups)

var.test(All_red46_mean, All_red7_mean)
t.test(All_red46_mean, All_red7_mean)

#Means testing for normality
shapiro.test(All_red46_mean)
shapiro.test(All_red7_mean)

skewness(All_red46_mean)
kurtosis(All_red46_mean)

skewness(All_red7_mean)
kurtosis(All_red7_mean)

shapiro.test(mean_of_all_groups)
wilcox.test(mean_of_all_groups)

shapiro.test(mean_of_complex_groups)
wilcox.test(mean_of_free_groups, mean_of_complex_groups, paired=TRUE)

shapiro.test(mean_of_free_groups)
wilcox.test(mean_of_free_groups)

t.test(mean_of_complex_groups, mean_of_free_groups)


#Only works if data is normally distributed but is under 30, hard to tell
t.test(range_of_free_groups, range_of_complex_groups, paired=TRUE)

show(range_10plus_redundgroups)
show(sd_10plus_redundgroups)
