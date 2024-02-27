load("data-danny-same-other-cv-01-19.RData")
library(data.table)
score.wide <- dcast(
  score.dt,
  train.groups + test.fold + test.group + task_id ~ algorithm,
  value.var="regr.mse")
score.wide[is.na(cv_glmnet), cv_glmnet := featureless]
score.tall <- melt(
  score.wide,
  measure=c("cv_glmnet","featureless"),
  variable.name="algorithm",
  value.name="regr.mse")


library(ggplot2)
gg <- ggplot()+
  geom_text(aes(
    regr.mse, train.groups, label=train.groups),
    vjust=-0.5,
    hjust=1,
    data=data.table(regr.mse=Inf, train.groups=c("same","other","all")))+
  geom_point(aes(
    regr.mse, train.groups, color=algorithm),
    shape=1,
    data=score.tall)+
  facet_grid(test.group ~ task_id, scales="free")+
  scale_x_log10()
png("data-danny-same-other-cv-figure-01-19.png", height=5, width=60, units="in", res=100)
print(gg)
dev.off()

## TODO compute p-value for difference between cv_glmnet and featureless
gg2 <- ggplot()+
  geom_point(aes(
    cv_glmnet-featureless, train.groups),
    data=score.wide)+
  facet_grid(. ~ test.group, labeller=label_both)