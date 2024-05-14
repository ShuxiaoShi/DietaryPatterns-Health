library(fda)
library(MASS)

fpca <- function(A,rank=0,ifscale=TRUE,lambda=NULL){
  if(ifscale){A <- scale(as.matrix(A))[,]}
  A.svd <- svd(A)
  if(rank==0){
    d <- A.svd$d
  } else {
    d <- A.svd$d-A.svd$d[min(rank+1,nrow(A),ncol(A))]
  }
  d <- d[d > 1e-8]
  r <- length(d)
  prop <- d^2; info <- sum(prop)/sum(A.svd$d^2);prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  rlt <- list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop,info=info)
  return(rlt)
}
data <- read.csv("food_class2.csv")

dt<-as.matrix(scale(data[,2:dim(data)[2]])) 
rlt<-fpca(dt,rank=3)
dt_new<-cbind(data[,1],rlt$X)
colnames(dt_new)<-c("n_eid","pc1","pc2","pc3")
dt_new_raw <- as.data.frame(dt_new)
write.csv(dt_new_raw,file="diet_pca_raw.csv",row.names = FALSE)


dt_new2 <- as.data.frame(dt_new)
summary(dt_new2$pc1)
sd(dt_new2$pc1)
dt_new2$pc1_scale <- scale(dt_new2$pc1,center = FALSE,scale = TRUE)
dt_new2$pc2_scale <- scale(dt_new2$pc2,center = FALSE,scale = TRUE)
dt_new2$pc3_scale <- scale(dt_new2$pc3,center = FALSE,scale = TRUE)

write.csv(dt_new2[-2:-4],file="diet_pca_white.csv",row.names = FALSE)

write.csv(dt_new2[-2:-4],file="C:/Users/15463/Desktop/BaiduSyncdisk/lefse/LEMMA/paper/journal/NF/first_response/data/diet_pca_white.csv",row.names = FALSE)
#dt_fifth$pc1_category <- cut(dt_fifth$pc1, breaks = quantile(dt_fifth$pc1, probs = 0:5/5), include.lowest = TRUE, labels = 1:5)
#data$pc1_category <- as.factor(data$pc1_category)
#dummy_vars <- model.matrix(~ pc1_category - 1, data = data)
#data <- cbind(data, dummy_vars)
dt_fifth <- as.data.frame(dt_new)
create_dummy_from_quantiles <- function(data, var_name, dummy_prefix) {
  # 检查变量名是否在数据框中
  if (!var_name %in% names(data)) {
    stop("Variable not found in the data frame")
  }
  
  # 按分位数划分变量
  quantile_breaks <- quantile(data[[var_name]], probs = 0:5/5, na.rm = TRUE)
  category_var <- cut(data[[var_name]], breaks = quantile_breaks, include.lowest = TRUE, labels = 1:5)
  
  # 转换为因子
  category_var <- as.factor(category_var)
  
  # 创建哑变量
  dummy_vars <- model.matrix(~ category_var - 1)
  colnames(dummy_vars) <- paste(dummy_prefix, colnames(dummy_vars), sep = "_")
  
  # 合并哑变量到原始数据框
  data <- cbind(data, dummy_vars)
  
  return(data)
}


dt_fifth_updated <- create_dummy_from_quantiles(dt_fifth, 'pc1', 'pc1')
dt_fifth_updated <- create_dummy_from_quantiles(dt_fifth_updated, 'pc2', 'pc2')
dt_fifth_updated <- create_dummy_from_quantiles(dt_fifth_updated, 'pc3', 'pc3')

dt_fifth_final <- dt_fifth_updated[-2:-4]
write.csv(dt_fifth_final,file="diet_pca_fifth.csv",row.names = FALSE)
###other methods

#######summary 0 value#####
library(dplyr)
library(purrr)
diet_pca <- read.csv("C:/Users/15463/Desktop/BaiduSyncdisk/lefse/LEMMA/pca_diet/food_class2.csv")
diet_pca <- read.csv("C:/Users/15463/Desktop/BaiduSyncdisk/lefse/LEMMA/food_group93/food93_group.csv")
# 除去第一列
diet_pca_no_id <- diet_pca[, -1]

# 计算每列中0的个数
zero_counts <- sapply(diet_pca_no_id, function(x) sum(x == 0, na.rm = TRUE))

# 计算每列中0的占比
zero_proportions <- sapply(diet_pca_no_id, function(x) mean(x == 0, na.rm = TRUE))

# 创建一个汇总的DataFrame
summary_df <- data.frame(
  Variable = names(zero_counts),
  Zero_Count = zero_counts,
  Zero_Proportion = zero_proportions
)


write.csv(summary_df,file="C:/Users/15463/Desktop/BaiduSyncdisk/lefse/LEMMA/paper/journal/NF/first_response/data/summary_zero.csv",row.names = FALSE)



