wilc_test <- wilcox.test(nonredund_complex_pangles, nonredund_free_pangles)
print(wilc_test)

kolgotest<-ks.test(nonredund_complex_pangles, nonredund_free_pangles)
print(kolgotest)

variation_test <- var.test(nonredund_complex_pangles, nonredund_free_pangles)
print(variation_test)



