# 清空环境
gc()
rm(list = ls())

# 加载所需包
library(dplyr)
library(tidyr)
library(preprocessCore)  # 用于分位数归一化
library(sva)           # 用于ComBat批次效应校正
library(limma)         # 用于线性回归归一化

# 设置工作目录
setwd("D:/BiologicalCodeinR/proteomeExpert-proteomeExpert/FFPE_glycopeptides/FFPE_glycopeptides")

# 加载数据
data2_1=read.csv("FFPE_mix2_glyco_4h_3ug_60K_repeat1.csv",header = T,sep = ",")
data2_1$raw.file=rep("raw2_1",as.numeric(length(data2_1$MS1)))

data2_2=read.csv("FFPE_mix2_glyco_4h_3ug_60K_repeat2.csv",header = T,sep = ",")
data2_2$raw.file=rep("raw2_2",as.numeric(length(data2_2$MS1)))

data2_3=read.csv("FFPE_mix2_glyco_4h_3ug_60K_repeat3.csv",header = T,sep = ",")
data2_3$raw.file=rep("raw2_3",as.numeric(length(data2_3$MS1)))

data3_1=read.csv("FFPE_mix3_glyco_4h_3ug_60K_repeat1.csv",header = T,sep = ",")
data3_1$raw.file=rep("raw3_1",as.numeric(length(data3_1$MS1)))

data3_2=read.csv("FFPE_mix3_glyco_4h_3ug_60K_repeat2.csv",header = T,sep = ",")
data3_2$raw.file=rep("raw3_2",as.numeric(length(data3_2$MS1)))

data3_3=read.csv("FFPE_mix3_glyco_4h_3ug_60K_repeat3.csv",header = T,sep = ",")
data3_3$raw.file=rep("raw3_3",as.numeric(length(data3_3$MS1)))

data4_1=read.csv("FFPE_mix4_glyco_4h_3ug_60K_repeat1.csv",header = T,sep = ",")
data4_1$raw.file=rep("raw4_1",as.numeric(length(data4_1$MS1)))

data4_2=read.csv("FFPE_mix4_glyco_4h_3ug_60K_repeat2.csv",header = T,sep = ",")
data4_2$raw.file=rep("raw4_2",as.numeric(length(data4_2$MS1)))

data4_3=read.csv("FFPE_mix4_glyco_4h_3ug_60K_repeat3.csv",header = T,sep = ",")
data4_3$raw.file=rep("raw4_3",as.numeric(length(data4_3$MS1)))

data5_1=read.csv("FFPE_mix5_glyco_4h_3ug_60K_repeat1.csv",header = T,sep = ",")
data5_1$raw.file=rep("raw5_1",as.numeric(length(data5_1$MS1)))

data5_2=read.csv("FFPE_mix5_glyco_4h_3ug_60K_repeat2.csv",header = T,sep = ",")
data5_2$raw.file=rep("raw5_2",as.numeric(length(data5_2$MS1)))

data5_3=read.csv("FFPE_mix5_glyco_4h_3ug_60K_repeat3.csv",header = T,sep = ",")
data5_3$raw.file=rep("raw5_3",as.numeric(length(data5_3$MS1)))

data6_1=read.csv("FFPE_mix6_glyco_4h_3ug_60K_repeat1.csv",header = T,sep = ",")
data6_1$raw.file=rep("raw6_1",as.numeric(length(data6_1$MS1)))

data6_2=read.csv("FFPE_mix6_glyco_4h_3ug_60K_repeat2.csv",header = T,sep = ",")
data6_2$raw.file=rep("raw6_2",as.numeric(length(data6_2$MS1)))

data6_3=read.csv("FFPE_mix6_glyco_4h_3ug_60K_repeat3.csv",header = T,sep = ",")
data6_3$raw.file=rep("raw6_3",as.numeric(length(data6_3$MS1)))

data7_1=read.csv("FFPE_mix7_glyco_4h_3ug_60K_repeat1.csv",header = T,sep = ",")
data7_1$raw.file=rep("raw7_1",as.numeric(length(data7_1$MS1)))

data7_2=read.csv("FFPE_mix7_glyco_4h_3ug_60K_repeat2.csv",header = T,sep = ",")
data7_2$raw.file=rep("raw7_2",as.numeric(length(data7_2$MS1)))

data7_3=read.csv("FFPE_mix7_glyco_4h_3ug_60K_repeat3.csv",header = T,sep = ",")
data7_3$raw.file=rep("raw7_3",as.numeric(length(data7_3$MS1)))

data8_1=read.csv("FFPE_mix8_glyco_4h_3ug_60K_repeat1.csv",header = T,sep = ",")
data8_1$raw.file=rep("raw8_1",as.numeric(length(data8_1$MS1)))

data8_2=read.csv("FFPE_mix8_glyco_4h_3ug_60K_repeat2.csv",header = T,sep = ",")
data8_2$raw.file=rep("raw8_2",as.numeric(length(data8_2$MS1)))

data8_3=read.csv("FFPE_mix8_glyco_4h_3ug_60K_repeat3.csv",header = T,sep = ",")
data8_3$raw.file=rep("raw8_3",as.numeric(length(data8_3$MS1)))

data9_1=read.csv("FFPE_mix9_glyco_4h_3ug_60K_repeat1.csv",header = T,sep = ",")
data9_1$raw.file=rep("raw9_1",as.numeric(length(data9_1$MS1)))

data9_2=read.csv("FFPE_mix9_glyco_4h_3ug_60K_repeat2.csv",header = T,sep = ",")
data9_2$raw.file=rep("raw9_2",as.numeric(length(data9_2$MS1)))

data9_3=read.csv("FFPE_mix9_glyco_4h_3ug_60K_repeat3.csv",header = T,sep = ",")
data9_3$raw.file=rep("raw9_3",as.numeric(length(data9_3$MS1)))


# 合并所有数据
data_all <- rbind(data2_1, data2_2, data2_3, data3_1, data3_2, data3_3,
                  data4_1, data4_2, data4_3, data5_1, data5_2, data5_3,
                  data6_1, data6_2, data6_3, data7_1, data7_2, data7_3,
                  data8_1, data8_2, data8_3, data9_1, data9_2, data9_3)

# 按原始文件拆分数据
list_all <- split(data_all, data_all$raw.file)


### 数据预处理函数
preprocess_data <- function(data) {
  # 初步筛选
  data <- data %>%
    filter(FDR <= 0.01) %>%
    filter(NGLYCAN != "NA") %>%
    filter(Protein.Name != "DECOY") %>%
    separate(Reporters, sep = ";", into = paste0("R_", 1:10)) %>%
    separate("R_1", sep = "/", into = c("R_126", "intensity_126")) %>%
    separate("R_2", sep = "/", into = c("R_127C", "intensity_127C")) %>%
    separate("R_3", sep = "/", into = c("R_127N", "intensity_127N")) %>%
    separate("R_4", sep = "/", into = c("R_128C", "intensity_128C")) %>%
    separate("R_5", sep = "/", into = c("R_128N", "intensity_128N")) %>%
    separate("R_6", sep = "/", into = c("R_129C", "intensity_129C")) %>%
    separate("R_7", sep = "/", into = c("R_129N", "intensity_129N")) %>%
    separate("R_8", sep = "/", into = c("R_130C", "intensity_130C")) %>%
    separate("R_9", sep = "/", into = c("R_130N", "intensity_130N")) %>%
    separate("R_10", sep = "/", into = c("R_131", "intensity_131")) %>%
    select(-starts_with("R_")) %>%
    select(Protein.Accession, Sequence, NGLYCAN, raw.file, starts_with("intensity_")) %>%
    separate(NGLYCAN, sep = ":", into = c("NGLYCAN", "NO")) %>%
    select(-NO) %>%
    unite(Sequence, NGLYCAN, col = "Glycopeptides", sep = "-") %>%
    mutate(Protein.Accession = sub(";.*", "", Protein.Accession))
  
  # 转换强度列为数值型
  exp_cols <- grep("^intensity", names(data), value = TRUE)
  data[exp_cols] <- lapply(data[exp_cols], as.numeric)
  
  # 过滤掉强度为0或NA的值
  data <- data %>% filter_at(vars(exp_cols), all_vars(. > 0))
  
  return(data)
}


### 组内归一化方法

# 1. 中位数归一化 (原代码中的方法)
median_normalization <- function(data) {
  # 转换为log2尺度
  exp_cols <- grep("^intensity", names(data), value = TRUE)
  data[exp_cols] <- lapply(data[exp_cols], function(x) log2(x + 1))
  
  # 长格式转换
  data_long <- gather(data, key = sample, value = intensity, -c(Protein.Accession, Glycopeptides, raw.file))
  
  # 对每个样本的每个糖肽计算中位数
  normalized_data <- data_long %>%
    group_by(sample, Glycopeptides) %>%
    mutate(normalized_intensity = intensity - median(intensity, na.rm = TRUE)) %>%
    ungroup()
  
  return(normalized_data)
}

# 2. 总离子强度归一化 (TIC Normalization)
tic_normalization <- function(data) {
  # 转换为log2尺度
  exp_cols <- grep("^intensity", names(data), value = TRUE)
  data[exp_cols] <- lapply(data[exp_cols], function(x) log2(x + 1))
  
  # 计算每个样本的TIC
  data$tic <- rowSums(data[exp_cols])
  
  # 计算所有样本的平均TIC
  mean_tic <- mean(data$tic)
  
  # 归一化每个强度值
  data[exp_cols] <- lapply(exp_cols, function(col) {
    data[[col]] * mean_tic / data$tic
  })
  
  # 长格式转换
  data_long <- gather(data, key = sample, value = intensity, -c(Protein.Accession, Glycopeptides, raw.file, tic))
  
  return(data_long)
}

# 3. 分位数归一化 (Quantile Normalization)
quantile_normalization <- function(data) {
  # 提取强度数据
  exp_cols <- grep("^intensity", names(data), value = TRUE)
  intensity_matrix <- as.matrix(data[exp_cols])
  
  # 分位数归一化
  normalized_matrix <- normalize.quantiles(intensity_matrix)
  
  # 转换回数据框
  normalized_data <- data.frame(data[, c("Protein.Accession", "Glycopeptides", "raw.file")])
  normalized_data[exp_cols] <- normalized_matrix
  
  # 转换为log2尺度
  normalized_data[exp_cols] <- lapply(normalized_data[exp_cols], function(x) log2(x + 1))
  
  # 长格式转换
  data_long <- gather(normalized_data, key = sample, value = intensity, -c(Protein.Accession, Glycopeptides, raw.file))
  
  return(data_long)
}

# 4. 概率商归一化 (PQN)
pqn_normalization <- function(data) {
  # 提取强度数据
  exp_cols <- grep("^intensity", names(data), value = TRUE)
  intensity_matrix <- as.matrix(data[exp_cols])
  
  # 转换为log2尺度
  intensity_matrix <- log2(intensity_matrix + 1)
  
  # 概率商归一化
  ref_sample <- rowMedians(intensity_matrix)
  pqn_matrix <- t(apply(intensity_matrix, 1, function(row) {
    quotients <- row / ref_sample
    median_quotient <- median(quotients, na.rm = TRUE)
    row - log2(median_quotient)
  }))
  
  # 转换回数据框
  normalized_data <- data.frame(data[, c("Protein.Accession", "Glycopeptides", "raw.file")])
  normalized_data[exp_cols] <- pqn_matrix
  
  # 长格式转换
  data_long <- gather(normalized_data, key = sample, value = intensity, -c(Protein.Accession, Glycopeptides, raw.file))
  
  return(data_long)
}


### 组间归一化方法

# 1. 看家蛋白归一化 (假设已知看家蛋白)
housekeeping_normalization <- function(data, hk_proteins = c("P54855", "P00738")) {
  # 筛选看家蛋白
  hk_data <- data %>% filter(Protein.Accession %in% hk_proteins)
  
  # 计算每个样本中看家蛋白的平均强度
  hk_means <- hk_data %>%
    group_by(raw.file, sample) %>%
    summarise(hk_mean = mean(intensity, na.rm = TRUE), .groups = 'drop') %>%
    ungroup()
  
  # 计算所有样本的平均看家蛋白强度
  overall_hk_mean <- mean(hk_means$hk_mean)
  
  # 计算校正因子
  hk_means$correction_factor <- overall_hk_mean / hk_means$hk_mean
  
  # 应用校正因子
  normalized_data <- data %>%
    left_join(hk_means, by = c("raw.file", "sample")) %>%
    mutate(normalized_intensity = intensity * correction_factor) %>%
    select(-hk_mean, -correction_factor)
  
  return(normalized_data)
}

# 2. 线性回归归一化
linear_regression_normalization <- function(data) {
  # 创建设计矩阵
  samples <- unique(data$sample)
  design <- model.matrix(~ 0 + samples)
  
  # 提取表达矩阵
  exp_matrix <- data %>%
    select(Glycopeptides, sample, intensity) %>%
    group_by(Glycopeptides, sample) %>%
    summarise(intensity = mean(intensity, na.rm = TRUE), .groups = 'drop') %>%
    ungroup() %>%
    spread(key = sample, value = intensity) %>%
    as.data.frame()
  
  rownames(exp_matrix) <- exp_matrix$Glycopeptides
  exp_matrix$Glycopeptides <- NULL
  
  # 进行线性回归拟合
  fit <- lmFit(as.matrix(exp_matrix), design)
  
  # 计算归一化后的表达矩阵
  normalized_matrix <- as.matrix(exp_matrix) - fit$coefficients %*% t(design)
  
  # 转换回长格式
  normalized_data <- as.data.frame(normalized_matrix) %>%
    tibble::rownames_to_column(var = "Glycopeptides") %>%
    tidyr::gather(key = "sample", value = "normalized_intensity", -Glycopeptides) %>%
    left_join(data %>% select(Glycopeptides, Protein.Accession, raw.file) %>% distinct(), by = "Glycopeptides")
  
  return(normalized_data)
}

# 3. ComBat批次效应校正
combat_normalization <- function(data) {
  # 创建批次向量 (假设raw.file对应批次)
  batches <- unique(data$raw.file)
  batch_vector <- as.numeric(factor(data$raw.file))
  
  # 创建表达矩阵
  exp_matrix <- data %>%
    select(Glycopeptides, sample, intensity) %>%
    group_by(Glycopeptides, sample) %>%
    summarise(intensity = mean(intensity, na.rm = TRUE), .groups = 'drop') %>%
    ungroup() %>%
    spread(key = sample, value = intensity) %>%
    as.data.frame()
  
  rownames(exp_matrix) <- exp_matrix$Glycopeptides
  exp_matrix$Glycopeptides <- NULL
  
  # 应用ComBat校正
  adjusted_matrix <- ComBat(as.matrix(exp_matrix), batch = batch_vector)
  
  # 转换回长格式
  normalized_data <- adjusted_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = "Glycopeptides") %>%
    gather(key = "sample", value = "normalized_intensity", -Glycopeptides) %>%
    left_join(data %>% select(Glycopeptides, Protein.Accession, raw.file) %>% distinct(), by = "Glycopeptides")
  
  return(normalized_data)
}


### 执行归一化流程

# ===== 归一化方法设置（如需修改归一化方法，请修改此处） =====
selected_zunei_method <- 'median'   # 组内归一化方法，可选：median, quantile, tic, pqn
selected_zujian_method <- 'channel_131'  # 组间归一化方法，可选：global_median, channel_131, linear
# ===== 归一化方法设置结束 =====

# 应用组内归一化
apply_zunei_normalization <- function(data, method) {
  preprocessed_data <- preprocess_data(data)
  
  switch(method,
         "median" = median_normalization(preprocessed_data),
         "tic" = tic_normalization(preprocessed_data),
         "quantile" = quantile_normalization(preprocessed_data),
         "pqn" = pqn_normalization(preprocessed_data),
         stop("Unsupported normalization method")
  )
}

# 应用组间归一化
apply_zujian_normalization <- function(data, method) {
  switch(method,
         "housekeeping" = housekeeping_normalization(data),
         "linear" = linear_regression_normalization(data),
         "combat" = combat_normalization(data),
         stop("Unsupported normalization method")
  )
}

# 归一化流程
normalized_data_list <- lapply(list_all, function(data) {
  zunei_data <- apply_zunei_normalization(data, selected_zunei_method)
  zujian_data <- apply_zujian_normalization(zunei_data, selected_zujian_method)
  return(zujian_data)
})

final_normalized_data <- bind_rows(normalized_data_list)

# 保存结果到设定目录
output_file <- file.path(getwd(), paste0("normalized_data_", selected_zunei_method, "_", selected_zujian_method, ".csv"))
write.csv(final_normalized_data, file = output_file, row.names = FALSE)
cat("归一化结果已保存到：", output_file, "\n")

