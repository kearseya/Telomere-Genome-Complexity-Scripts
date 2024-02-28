library(tidyverse)
library(gridExtra)

l <- c('11_0.gc', '11_0.kc_perc', '13_0.fragments', '13_0.kc_len_cov',
       '13_0.new', '13_0.num_reads', '13_0.scqc', '2_0.scqc', '44_0.adj_cov',
       '44_0.kc_len_reg', '44_0.scqc', '45_0.avg_read_len', '45_0.gc',
       '45_0.kc_len', '45_0.len_reg', '4_0.gc', '4_0.kc_perc', '9_0.new',
       'f_0.kc_len_cov', 'f_0.kc_len_reg', 'f_0.new', 'len_cov', 'r_0.adj_cov',
       'short.len_reg', 'short.length', 'short.mapq', 'short.num_regs')

train <- read_csv("/path/to/project/TL_prediction/predictions_for_r_updated.csv")
train

test <- read_csv("/path/to/project/TL_prediction/GC_cov_testing/lets_see.csv")
test


transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column(.data = .) %>%
    tibble::as_tibble(.)
  return(t_df)
}


train_cov <- read_csv("/path/to/project/TL_prediction/coverages/coverage_for_r.csv")
train_cov <- transpose_df(train_cov)
test_cov <- read_csv("/path/to/project/TL_prediction/GC_cov_testing/coverages.csv")
test_cov <- transpose_df(test_cov)

train_cov <- train_cov %>% rename("sample"="rowname", "cov"="1")
test_cov <- test_cov %>% rename("sample"="rowname", "cov"="1")

train <- left_join(train, train_cov)
test <- left_join(test, test_cov)

train <- train %>% add_column(train_data=1)
test <- test %>% add_column(train_data=0)

all <- full_join(train, test)

all <- rowid_to_column(all, "X1")
all %<% m

plot(test$`13.num_reads`*test$cov)
plot(train$`13.num_reads`*train$cov)

plot(train$`13_0.num_reads`)
plot(test$`13_0.num_reads`)

##
bgi <- all %>% filter(train_data == 1) %>% select(total_reads)
ill <- all %>% filter(train_data == 0) %>% select(total_reads)

ggplot(all, aes(x=total_reads)) +
  geom_histogram(data=subset(all,train_data == 0),fill = "red", alpha = 0.5) +
  geom_histogram(data=subset(all,train_data == 1),fill = "blue", alpha = 0.5)

plot(all$cov, all$f_0.adj_cov*all$cov)




all$total_new

ggplot(all, aes(X1, `f_0.num_reads`, col=train_data)) + 
  geom_point() +
  theme(legend.position="none")
ggplot(all, aes(X1, `r_0.num_reads`, col=train_data)) + 
  geom_point() +
  theme(legend.position="none")
ggplot(all, aes(X1, `11_0.num_reads`, col=train_data)) + 
  geom_point() +
  theme(legend.position="none")

ggplot(all, aes(X1, `total_new`, col=train_data)) + 
  geom_point() +
  theme(legend.position="none")


ggplot(all, aes(X1, `13_0.num_reads`, col=train_data)) + 
  geom_point() +
  theme(legend.position="none")
ggplot(all, aes(X1, `44_0.num_reads`, col=train_data)) + 
  geom_point() +
  theme(legend.position="none")
ggplot(all, aes(X1, `45_0.num_reads`, col=train_data)) + 
  geom_point() +
  theme(legend.position="none")

ggplot(all, aes(X1, `mapq0.num_reads`, col=train_data)) + 
  geom_point() +
  theme(legend.position="none")
ggplot(all, aes(X1, `mapq1.num_reads`, col=train_data)) + 
  geom_point() +
  theme(legend.position="none")
ggplot(all, aes(X1, `total_reads`, col=train_data)) + 
  geom_point() +
  theme(legend.position="none")

all <- all %>% mutate(`total_ex_reads`=(`total_reads`-`mapq0.num_reads`-`mapq1.num_reads`-`short.num_reads`))

ggplot(all, aes(cov, (`total_ex_reads`*cov)/`total_adj_cov`, col=train_data)) + 
  geom_point() +
  geom_smooth() +
  ylab("Total raw reads / total adj cov") +
  theme(legend.position="none")

ggplot(all, aes(cov, `f_0.adj_cov`*cov, col=train_data)) + 
  geom_point() +
  #geom_smooth() +
  #ylab("Total raw reads / total adj cov") +
  theme(legend.position="none")

ggplot(all, aes(cov, `total_ex_reads`, col=train_data)) + 
  geom_point() +
  #geom_smooth() +
  theme(legend.position="none")


ggplot(all, aes(cov, `total_adj_cov`, col=train_data)) + 
  geom_point() +
  geom_smooth() +
  theme(legend.position="none")

ggplot(all, aes(cov, `total_kc_length`, col=train_data)) + 
  geom_point() +
  geom_smooth() +
  theme(legend.position="none")

ggplot(all, aes(cov, `total_length`, col=train_data)) + 
  geom_point() +
  geom_smooth() +
  theme(legend.position="none")

ggplot(all, aes(cov, `f_0.adj_cov`, col=train_data)) + 
  geom_point() +
  theme(legend.position="none")

ggplot(all, aes(cov, `r_0.adj_cov`, col=train_data)) + 
  geom_point() +
  theme(legend.position="none")

ggplot(all, aes(cov, `r_0.adj_cov`, col=train_data)) + 
  geom_point() +
  theme(legend.position="none")

ggplot(all, aes(cov, `44_0.adj_cov`, col=train_data)) + 
  geom_point() +
  theme(legend.position="none")

####### ADJUSTED COVERAGE
## same
c(13, 2, ) 
## intersect
c(45, 44)
## diff
c(9, 4) 
## curve
c(11)

ggplot(all, aes(cov, `f_0.kc_len`, col=train_data)) + 
  geom_point() +
  theme(legend.position="none")


plot(test$cov)


ggplot(all, aes(cov, `r_0.num_reads`/`r_0.adj_cov`, col=train_data)) + 
  geom_point(aes(col=short)) +
  theme(legend.position="none")


z_scores_r <- (train$total_reads-mean(train$total_reads))/sd(train$total_reads)
z_scores_e <- (test$total_reads-mean(test$total_reads))/sd(test$total_reads)

cv_scores_r <- (sd(train$total_reads)/mean(train$total_reads))
cv_scores_e <- (sd(test$total_reads)/mean(test$total_reads))

hist(cv_scores_r)
hist(cv_scores_e)
  
cv_scores_r
cv_scores_e

mean(z_scores_e)

all$short

ggplot(all, aes(cov, `13_0.num_reads`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)
ggplot(all, aes(cov, `2_0.scqc`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)

## messy upward
ggplot(all, aes(cov, `44_0.adj_cov`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)
## clean downward
ggplot(all, aes(cov, `f_0.adj_cov`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)

ggplot(all, aes(cov, `len_cov`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)
ggplot(all, aes(cov, `4_0.gc`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)
ggplot(all, aes(cov, `45_0.kc_len`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)

ggplot(all, aes(cov, `45_0.avg_read_len`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)


ggplot(all, aes(cov, `4_0.kc_perc`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)

ggplot(all, aes(cov, `f_0.kc_len_cov`)) + geom_point() + geom_smooth(method="lm") + geom_vline(xintercept=21)


ggplot(all, aes(cov, `f_0.gc`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)
ggplot(all, aes(cov, `f_0.adj_cov`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)
ggplot(all, aes(cov, `f_0.kc_len`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)

ggplot(all, aes(cov, `45_0.avg_read_len`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)



## correlated down curve
ggplot(all, aes(cov, `11_0.adj_cov`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)
ggplot(all, aes(cov, `f_0.adj_cov`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)
ggplot(all, aes(cov, `r_0.adj_cov`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)
## messy up curve
ggplot(all, aes(cov, `9_0.adj_cov`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)
ggplot(all, aes(cov, `45_0.adj_cov`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)
ggplot(all, aes(cov, `44_0.adj_cov`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)
ggplot(all, aes(cov, `2_0.adj_cov`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)
ggplot(all, aes(cov, `13_0.adj_cov`)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)


ggplot(train, aes(x=short, y=`11_0.gc`)) + 
  geom_boxplot()

regcov <- all$cov*all$f_0.adj_cov
plot(all$cov, regcov)

## all$forwardcov <- 
all <- mutate(all, forwardcov = cov*`f_0.adj_cov`) 
ggplot(all, aes(cov, forwardcov)) + geom_point(aes(col=short)) + geom_smooth(method="lm") + geom_vline(xintercept=21)

for (i in l){
  print(i)
  p <- ggplot(train, aes_string(x="short", y=train$i)) + 
    geom_boxplot()
}
p

p <- ggplot(train, aes_string(x="short", y=train$"11_0.gc")) + 
  geom_boxplot()
p

i <- 0
train_plot_list = list()
for (c in l) {
  print(c)
  p <- ggplot(train, aes_string(x="short", y=str(c))) + 
    geom_boxplot()
  train_plot_list[[i]] <- p
  i <- i+1
}

train_plot_list[[0]]


train_short <- train %>% filter(short==TRUE)
train_long <- train %>% filter(short==FALSE)

i <- 0
train_plot_list = list()
for (c in l) {
  print(c)
  p <- boxplot(train_long[c])
  train_plot_list[[i]] <- p
  i <- i+1
}

train_plot_list[[4]]


?boxplot




r1 <- ggplot(train, aes(x=short, y=`11_0.gc`)) + 
  geom_boxplot()
r2 <- ggplot(train, aes(x=short, y=`11_0.kc_perc`)) + 
  geom_boxplot()
r3 <- ggplot(train, aes(x=short, y=`13_0.fragments`)) + 
  geom_boxplot()
r4 <- ggplot(train, aes(x=short, y=`13_0.kc_len_cov`)) + 
  geom_boxplot() +
  ylim(0, 40000)
r5 <- ggplot(train, aes(x=short, y=`13_0.new`)) + 
  geom_boxplot()
r6 <- ggplot(train, aes(x=short, y=`13_0.num_reads`)) + 
  geom_boxplot()
r7 <- ggplot(train, aes(x=short, y=`13_0.scqc`)) + 
  geom_boxplot()
r8 <- ggplot(train, aes(x=short, y=`2_0.scqc`)) + 
  geom_boxplot()
r9 <- ggplot(train, aes(x=short, y=`44_0.adj_cov`)) + 
  geom_boxplot()
r10 <- ggplot(train, aes(x=short, y=`44_0.kc_len_reg`)) + 
  geom_boxplot()
r11 <- ggplot(train, aes(x=short, y=`44_0.scqc`)) + 
  geom_boxplot()
r12 <- ggplot(train, aes(x=short, y=`45_0.avg_read_len`)) + 
  geom_boxplot()
r13 <- ggplot(train, aes(x=short, y=`45_0.gc`)) + 
  geom_boxplot()
r14 <- ggplot(train, aes(x=short, y=`45_0.kc_len`)) + 
  geom_boxplot()
r15 <- ggplot(train, aes(x=short, y=`45_0.len_reg`)) + 
  geom_boxplot()
r16 <- ggplot(train, aes(x=short, y=`4_0.gc`)) + 
  geom_boxplot()
r17 <- ggplot(train, aes(x=short, y=`4_0.kc_perc`)) + 
  geom_boxplot()
r18 <- ggplot(train, aes(x=short, y=`9_0.new`)) + 
  geom_boxplot()
r19 <- ggplot(train, aes(x=short, y=`f_0.kc_len_cov`)) + 
  geom_boxplot()
r20 <- ggplot(train, aes(x=short, y=`f_0.kc_len_reg`)) + 
  geom_boxplot()
r21 <- ggplot(train, aes(x=short, y=`f_0.new`)) + 
  geom_boxplot()
r22 <- ggplot(train, aes(x=short, y=`len_cov`)) + 
  geom_boxplot()
r23 <- ggplot(train, aes(x=short, y=`r_0.adj_cov`)) + 
  geom_boxplot()
r24 <- ggplot(train, aes(x=short, y=`short.len_reg`)) + 
  geom_boxplot()
r25 <- ggplot(train, aes(x=short, y=`short.length`)) + 
  geom_boxplot()
r26 <- ggplot(train, aes(x=short, y=`short.mapq`)) + 
  geom_boxplot()
r27 <- ggplot(train, aes(x=short, y=`short.num_regs`)) + 
  geom_boxplot()
ggplot(train, aes(x=short, y=``)) + 
  geom_boxplot()



e1 <- ggplot(test, aes(x=short, y=`11_0.gc`)) + 
  geom_boxplot()
e2 <- ggplot(test, aes(x=short, y=`11_0.kc_perc`)) + 
  geom_boxplot()
e3 <- ggplot(test, aes(x=short, y=`13_0.fragments`)) + 
  geom_boxplot()
e4 <- ggplot(test, aes(x=short, y=`13_0.kc_len_cov`)) + 
  geom_boxplot()
e5 <- ggplot(test, aes(x=short, y=`13_0.new`)) + 
  geom_boxplot()
e6 <- ggplot(test, aes(x=short, y=`13_0.num_reads`)) + 
  geom_boxplot()
e7 <- ggplot(test, aes(x=short, y=`13_0.scqc`)) + 
  geom_boxplot()
e8 <- ggplot(test, aes(x=short, y=`2_0.scqc`)) + 
  geom_boxplot()
e9 <- ggplot(test, aes(x=short, y=`44_0.adj_cov`)) + 
  geom_boxplot()
e10 <- ggplot(test, aes(x=short, y=`44_0.kc_len_reg`)) + 
  geom_boxplot()
e11 <- ggplot(test, aes(x=short, y=`44_0.scqc`)) + 
  geom_boxplot()
e12 <- ggplot(test, aes(x=short, y=`45_0.avg_read_len`)) + 
  geom_boxplot()
e13 <- ggplot(test, aes(x=short, y=`45_0.gc`)) + 
  geom_boxplot()
e14 <- ggplot(test, aes(x=short, y=`45_0.kc_len`)) + 
  geom_boxplot()
e15 <- ggplot(test, aes(x=short, y=`45_0.len_reg`)) + 
  geom_boxplot()
e16 <- ggplot(test, aes(x=short, y=`4_0.gc`)) + 
  geom_boxplot()
e17 <- ggplot(test, aes(x=short, y=`4_0.kc_perc`)) + 
  geom_boxplot()
e18 <- ggplot(test, aes(x=short, y=`9_0.new`)) + 
  geom_boxplot()
e19 <- ggplot(test, aes(x=short, y=`f_0.kc_len_cov`)) + 
  geom_boxplot()
e20 <- ggplot(test, aes(x=short, y=`f_0.kc_len_reg`)) + 
  geom_boxplot()
e21 <- ggplot(test, aes(x=short, y=`f_0.new`)) + 
  geom_boxplot()
e22 <- ggplot(test, aes(x=short, y=`len_cov`)) + 
  geom_boxplot()
e23 <- ggplot(test, aes(x=short, y=`r_0.adj_cov`)) + 
  geom_boxplot()
e24 <- ggplot(test, aes(x=short, y=`short.len_reg`)) + 
  geom_boxplot()
e25 <- ggplot(test, aes(x=short, y=`short.length`)) + 
  geom_boxplot()
e26 <- ggplot(test, aes(x=short, y=`short.mapq`)) + 
  geom_boxplot()
e27 <- ggplot(test, aes(x=short, y=`short.num_regs`)) + 
  geom_boxplot()



###
grid.arrange(r9, e9, ncol =2)



## only num read, kc len are higher
## reverse adj_cov lower, 44 adj_cov higher

1'11_0.gc', 
2'11_0.kc_perc', 
3'13_0.fragments', 
4'13_0.kc_len_cov',
5'13_0.new', 
6'13_0.num_reads', 
7'13_0.scqc', 
8'2_0.scqc', 
9'44_0.adj_cov',
10'44_0.kc_len_reg', 
11'44_0.scqc', 
12'45_0.avg_read_len', 
13'45_0.gc',
14'45_0.kc_len', 
15'45_0.len_reg', 
16'4_0.gc', 
17'4_0.kc_perc', 
18'9_0.new',
19'f_0.kc_len_cov', 
20'f_0.kc_len_reg', 
21'f_0.new', 
22'len_cov', 
23'r_0.adj_cov',
24'short.len_reg', 
25'short.length', 
26'short.mapq', 
27'short.num_regs'




