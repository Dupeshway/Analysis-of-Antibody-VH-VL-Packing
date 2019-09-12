#Dataframe analysis 10 plus redundant structures:

str(pangle.data)

abgroups_10plus_redundancies<-c(7, 8, 12, 13,
                                15, 20, 33, 40, 41,
                                51, 55, 59, 65)

spaced_10plus_redundancies<-c(NaN, NaN, NaN, NaN, NaN, NaN, 7, 8, NaN, 
                              NaN, NaN, 12, 13, NaN, 15,NaN, NaN, NaN, NaN,
                              20,NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 
                              NaN, NaN, NaN, NaN, 33,NaN, NaN, NaN, NaN,
                              NaN, NaN, 40, 41, NaN, NaN, NaN, NaN, NaN, NaN,
                              NaN, NaN, NaN, 51,NaN, NaN, NaN, 55,NaN, NaN,
                              NaN, 59,NaN, NaN, NaN, NaN, NaN, 65, NaN)

#ab groups only
mean_10plus_redundgroups<-mean_of_all_groups[abgroups_10plus_redundancies]
range_10plus_redundgroups<-range_of_all_groups[abgroups_10plus_redundancies]
sd_10plus_redundgroups<-sd_of_all_groups[abgroups_10plus_redundancies]

#spaced_10plus_redundancies
mean_spaced_10plus_redundgroups<-mean_of_all_groups[spaced_10plus_redundancies]
range_spaced_10plus_redundgroups<-range_of_all_groups[spaced_10plus_redundancies]
sd_spaced_10plus_redundgroups<-sd_of_all_groups[spaced_10plus_redundancies]

#10_redundances
range_redund_10plus <-dataframe_all_groups$range_redund_10.
sd_redund_10plus <- dataframe_all_groups$sd_redund_10.

show(mean_10plus_redundgroups)
str(mean_10plus_redundgroups)
str(range_10plus_redundgroups)
str(sd_10plus_redundgroups)

show(sd_of_all_groups)
show(sd_10plus_redundgroups)

show(range_of_all_groups)
show(range_10plus_redundgroups)

