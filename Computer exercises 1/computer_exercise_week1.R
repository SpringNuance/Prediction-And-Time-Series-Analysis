
library(ggplot2)


tobacco_data<- read.table('tobacco.txt', header=T, sep = '\t')

#a)
lm_ill <- lm(ILL ~ CONSUMPTION, data = tobacco_data)
summary(lm_ill)


