# Helper functions for FFPE glycopeptide analysis (Simplified version - Starting directly from inter-group normalized data)

# Calculate glycan molar mass
calculate_molar_mass <- function(glyco_str) {
  matches <- regmatches(glyco_str, regexec("N(\\d+)H(\\d+)F(\\d+)S(\\d+)G(\\d+)", glyco_str))
  if (length(matches[[1]]) < 6) return(NA)
  
  N <- as.integer(matches[[1]][2])
  H <- as.integer(matches[[1]][3])
  F <- as.integer(matches[[1]][4])
  S <- as.integer(matches[[1]][5])
  G <- as.integer(matches[[1]][6])
  
  molar_mass_GlcNAc <- 221.24
  molar_mass_Mannose <- 180.16
  molar_mass_Fucose <- 164.16
  molar_mass_Sialic_Acid <- 309.24
  molar_mass_Glucose <- 180.16
  
  total_mass <- (N * molar_mass_GlcNAc) + (H * molar_mass_Mannose) + 
    (F * molar_mass_Fucose) + (S * molar_mass_Sialic_Acid) + 
    (G * molar_mass_Glucose)
  
  glycosidic_bonds <- (N + H + F + S + G - 1)
  water_loss <- glycosidic_bonds * 18.015
  final_mass <- total_mass - water_loss
  
  return(final_mass)
}

# Define t-test function
perform_t_test <- function(df1, df2, glycan, glyco_mass, comparison) {
  # Only perform t-test when both groups have at least 2 non-NA observations
  if (sum(!is.na(df1$Value)) >= 2 && sum(!is.na(df2$Value)) >= 2) {
    t_test_result <- t.test(df1$Value, df2$Value)
    pval <- t_test_result$p.value
  } else {
    pval <- NA
  }
  data.frame(
    Glycan = glycan,
    Glyco_Molar_Mass = glyco_mass,
    p_value = pval,
    comparison = comparison
  )
}

# Glycomics analysis function (strictly replicating 20250521.R)
glycomics_analysis_20250521 <- function(zujian_data, patient_info) {
  # 还原强度
  zujian_data$final_scale_intensity <- 2^zujian_data$final_scale_intensity
  # 去除无关列
  data_for_glycomics <- zujian_data %>% dplyr::select(-Protein.Accession, -raw.file, -sample)
  # 拆分Glycopeptides
  data_for_glycomics <- data_for_glycomics %>%
    tidyr::separate(Glycopeptides, sep = "-", into = c("Peptides", "Glycan"), remove = FALSE) %>%
    dplyr::select(-Peptides)
  # 聚合
  data_for_glycomics <- data_for_glycomics %>%
    dplyr::group_by(Glycan, sample_id) %>%
    dplyr::summarise(total_intensity = sum(final_scale_intensity), .groups = 'drop')
  # 合并病人信息
  data_for_glycomics <- merge(data_for_glycomics, patient_info, by.x = "sample_id", by.y = "Sample")
  # 保留变量名一致
  data_for_glycomics <- dplyr::select(data_for_glycomics, -sample_id, -patients, -repeat.)
  # 宽表
  df_wide_for_glycomics <- data_for_glycomics %>%
    tidyr::pivot_wider(names_from = patients_repeat, values_from = total_intensity)
  # NA处理
  df_wide_for_glycomics[df_wide_for_glycomics == "NA"] <- NA
  # 有效值比例
  df_wide_for_glycomics <- df_wide_for_glycomics %>%
    dplyr::rowwise() %>%
    dplyr::mutate(valid_ratio = mean(as.numeric(!is.na(c_across(-Glycan)))))
  # 保留有效值大于0.5
  df_filtered <- df_wide_for_glycomics %>%
    dplyr::filter(valid_ratio > 0.5) %>%
    dplyr::select(-valid_ratio)
  # KNN填补
  df_filled <- VIM::kNN(df_filtered, variable = setdiff(colnames(df_filtered), "Glycan")) %>%
    dplyr::select(Glycan, dplyr::everything()) %>%
    dplyr::select(-matches("_imp$"))
  # 长表
  df_long_filled_glycomics <- df_filled %>%
    tidyr::pivot_longer(cols = -Glycan, names_to = "Sample", values_to = "Value")
  # 分割出cancertype和repeat，严格与20250521.R一致
  df_long_filled_glycomics <- df_long_filled_glycomics %>%
    tidyr::separate(Sample, sep = "\\.", into = c("patients", "cancertype")) %>%
    tidyr::separate(cancertype, sep = "_", into = c("cancertype", "repeat"))

  # 动态提取所有实际存在的cancertype
  cancertype_levels <- unique(df_long_filled_glycomics$cancertype)

  # 计算摩尔质量
  df_long_filled_glycomics <- df_long_filled_glycomics %>%
    dplyr::mutate(Glyco_Molar_Mass = sapply(Glycan, calculate_molar_mass))
  # 显著性差异
  # 动态分组
  cancertype_list <- split(df_long_filled_glycomics, df_long_filled_glycomics$cancertype)
  # 只做两两比较（如有多组，自动两两比较）
  cancertype_names <- names(cancertype_list)
  ttest_results <- list()
  if (length(cancertype_names) >= 2) {
    for (i in 1:(length(cancertype_names)-1)) {
      for (j in (i+1):length(cancertype_names)) {
        group1 <- cancertype_names[i]
        group2 <- cancertype_names[j]
        df1 <- cancertype_list[[group1]]
        df2 <- cancertype_list[[group2]]
        ttest_results[[paste(group1, group2, sep = " vs ")]] <- lapply(unique(df_long_filled_glycomics$Glycan), function(glycan) {
    glyco_mass <- unique(df_long_filled_glycomics$Glyco_Molar_Mass[df_long_filled_glycomics$Glycan == glycan])
          perform_t_test(df1[df1$Glycan == glycan, ], df2[df2$Glycan == glycan, ], glycan, glyco_mass, paste(group1, "vs", group2))
        })
      }
    }
  }
  results_df <- dplyr::bind_rows(unlist(ttest_results, recursive = FALSE))
  significant_results <- results_df %>% dplyr::filter(p_value < 0.05)
  # 归一化
  df_long_filled_glycomics_summary <- df_long_filled_glycomics %>%
    dplyr::group_by(Glyco_Molar_Mass, cancertype) %>%
    dplyr::summarise(mean_intensity = mean(Value, na.rm = TRUE), .groups = 'drop')
  data_for_normalized_glycomics <- df_long_filled_glycomics_summary %>%
    dplyr::group_by(cancertype) %>%
    dplyr::mutate(Mean_Value_normalized = (mean_intensity / max(mean_intensity)) * 100) %>%
    dplyr::ungroup()
  list(
    wide = df_wide_for_glycomics,
    long = df_long_filled_glycomics,
    normalized = data_for_normalized_glycomics,
    significant = significant_results,
    cancertype_levels = cancertype_levels
  )
}

# 糖肽定性分析函数（严格复刻20250521.R）
glycopeptide_analysis_20250521 <- function(zujian_data, patient_info, glycan_type) {
  zujian_data$final_scale_intensity <- 2^zujian_data$final_scale_intensity
  # 去除无关列
  data <- zujian_data %>% dplyr::select(-raw.file, -sample)
  # 宽表
  df_wide <- data %>% tidyr::pivot_wider(names_from = sample_id, values_from = final_scale_intensity)
  df_wide <- df_wide[, !grepl("131", names(df_wide))]
  df_wide$Glycopeptides <- paste(df_wide$Glycopeptides, df_wide$Protein.Accession, sep = "-")
  df_wide <- df_wide %>% dplyr::select(-Protein.Accession)
  # 列名替换为patients_repeat
  old_colnames <- colnames(df_wide)
  samples <- patient_info$Sample
  new_colnames <- patient_info$patients_repeat
  name_map <- setNames(new_colnames, samples)
  for (i in 1:length(old_colnames)) {
    if (old_colnames[i] %in% names(name_map)) {
      old_colnames[i] <- name_map[[old_colnames[i]]]
    }
  }
  colnames(df_wide) <- old_colnames
  # gather为长表
  df_long <- df_wide %>%
    tidyr::pivot_longer(
      cols = -c(Glycopeptides),
      names_to = "Sample",
      values_to = "Value"
    ) %>%
    dplyr::filter(!is.na(Value))
  # 合并糖型信息
  df_for_dingxing_analysis <- df_wide
  df_for_dingxing_analysis$glycopeptides <- df_for_dingxing_analysis$Glycopeptides
  df_for_dingxing_analysis <- df_for_dingxing_analysis %>%
    tidyr::separate(Glycopeptides, into = c("peptides", "glycan", "Protein"), sep = "-", extra = "merge")
  df_for_dingxing_analysis <- merge(df_for_dingxing_analysis, glycan_type)
  last_col <- tail(names(df_for_dingxing_analysis), n = 1)
  first_three_cols <- c("glycan", "peptides", "Protein")
  if (ncol(df_for_dingxing_analysis) > 4) {
    middle_cols <- names(df_for_dingxing_analysis)[4:(ncol(df_for_dingxing_analysis)-1)]
  } else {
    middle_cols <- character(0)
  }
  new_col_order <- c(first_three_cols, last_col, middle_cols)
  df_for_dingxing_analysis <- df_for_dingxing_analysis[, new_col_order]
  df_for_dingxing_analysis_long <- tidyr::gather(df_for_dingxing_analysis, key = "patient_info", value = "value", -glycan, -peptides, -Protein, -glycan.type,-glycopeptides)
  df_for_dingxing_analysis_long <- subset(df_for_dingxing_analysis_long, !is.na(df_for_dingxing_analysis_long$value))
  df_for_dingxing_analysis_long <- df_for_dingxing_analysis_long %>%
    tidyr::separate(patient_info, into = c("patients", "repeat"), sep = "_", extra = "merge") %>%
    tidyr::separate(patients, into = c("patients", "cancer_type"), sep = "\\.", extra = "merge")

  # 动态提取所有实际存在的cancer_type
  cancer_type_levels <- unique(df_for_dingxing_analysis_long$cancer_type)

  # unique统计
  unique_counts_by_cancer_type <- df_for_dingxing_analysis_long %>%
    dplyr::filter(value > 0) %>%
    dplyr::group_by(cancer_type) %>%
    dplyr::summarise(
      unique_proteins = n_distinct(Protein),
      unique_glycans = n_distinct(glycan),
      unique_peptides = n_distinct(peptides),
      unique_glycopeptides = n_distinct(glycopeptides)
    ) %>%
    tidyr::pivot_longer(cols = c(unique_proteins, unique_glycans, unique_peptides, unique_glycopeptides),
                        names_to = "unique_type", values_to = "count")
  peptide_glycan <- df_for_dingxing_analysis_long[, c("glycan", "peptides", "glycan.type")] %>% unique()
  glycan_counts_per_peptide <- peptide_glycan %>%
    dplyr::group_by(peptides) %>%
    dplyr::summarise(unique_glycans = n_distinct(glycan))
  glycopeptide_glycan <- peptide_glycan
  glycopeptide_glycan$glycopep <- paste(glycopeptide_glycan$peptides, glycopeptide_glycan$glycan, sep = "_")
  glycopeptide_glycan <- unique(glycopeptide_glycan)
  glycan_type_counts <- glycopeptide_glycan %>%
    dplyr::group_by(glycan.type) %>%
    dplyr::summarise(glycopeptides_count = n())
  glyco_glycan_sites <- merge(glycopeptide_glycan, as.data.frame(glycan_counts_per_peptide))
  colnames(glyco_glycan_sites) <- c("glycan", "peptides", "glycan.type", "glycopep", "sites")
  glyco_glycan_sites <- glyco_glycan_sites %>% dplyr::select(-glycopep)
  colnames(glyco_glycan_sites) <- c("peptides", "glycan", "glycan.type", "sites")
  
  # 新增：生成三个中间结果数据表
  # 1. 50%有效值保留_knn缺失值填补之前
  df_wide_for_knn <- df_wide
  df_wide_for_knn[df_wide_for_knn == "NA"] <- NA
  df_wide_for_knn <- df_wide_for_knn %>%
    dplyr::rowwise() %>%
    dplyr::mutate(valid_ratio = mean(as.numeric(!is.na(c_across(-contains("Glycopeptides"))))))
  df_filtered <- df_wide_for_knn %>%
    dplyr::filter(valid_ratio > 0.5) %>%
    dplyr::select(-valid_ratio)
  
  # 2. KNN填补
  df_filled <- VIM::kNN(df_filtered, variable = setdiff(colnames(df_filtered), "Glycopeptides"))
  df_filled <- dplyr::select(df_filled, -ends_with("_imp"))
  # 强制整理列名格式，保证与20250521.R一致
  colnames(df_filled) <- gsub("\\.+", "_", colnames(df_filled))
  
  # 3. 样品重复相关性
  data_for_cor <- df_filled
  rownames(data_for_cor) <- data_for_cor$Glycopeptides
  data_for_cor <- data_for_cor %>% dplyr::select(-Glycopeptides)
  data_for_cor <- data_for_cor[, order(names(data_for_cor))]
  correlation_matrix <- cor(data_for_cor, use = "pairwise.complete.obs")
  
  # 4. 填补缺失值_zscore后
  # 1. pivot_longer
  df_long_mean <- df_filled %>%
    tidyr::pivot_longer(
      cols = -Glycopeptides,
      names_to = c("Sample", "Repeat"),
      names_pattern = "(.*)_(repeat[0-9]+)"
    )
  # 2. group_by取均值
  df_mean <- df_long_mean %>%
    dplyr::group_by(Glycopeptides, Sample) %>%
    dplyr::summarize(Mean_Expression = mean(value, na.rm = TRUE), .groups = 'drop')
  # 3. pivot_wider
  df_wide_mean <- df_mean %>%
    tidyr::pivot_wider(
      names_from = Sample,
      values_from = Mean_Expression
    )
  # 4. 保证排序一致
  df_wide_mean <- df_wide_mean[order(df_wide_mean$Glycopeptides), ]
  sample_cols <- setdiff(colnames(df_wide_mean), "Glycopeptides")
  df_wide_mean <- df_wide_mean[, c("Glycopeptides", sort(sample_cols))]
  # 5. min-max归一化
  df_zscore <- df_wide_mean
  df_zscore[, -1] <- apply(df_wide_mean[, -1], 2, function(x) (x - min(x, na.rm=TRUE))/(max(x, na.rm=TRUE) - min(x, na.rm=TRUE)))
  
  list(
    wide = df_wide,
    long = df_long,
    unique_counts = unique_counts_by_cancer_type,
    peptide_glycan = peptide_glycan,
    glycan_counts_per_peptide = glycan_counts_per_peptide,
    glycan_type_counts = glycan_type_counts,
    glyco_glycan_sites = glyco_glycan_sites,
    df_for_dingxing_analysis_long = df_for_dingxing_analysis_long,
    knn_input = df_filtered,                # 50%有效值保留_knn缺失值填补之前
    correlation_matrix = correlation_matrix, # 样品重复相关性
    zscore = df_zscore,                      # 填补缺失值_zscore后
    cancer_type_levels = cancer_type_levels
  )
}

# 糖肽定量分析函数（缺失值处理）
handle_missing_values_glycopeptide <- function(data, patient_info, valid_ratio = 0.1, use_knn = TRUE, knn_k = 5) {
  # 将数据框转换为宽格式
  df_wide <- data %>%
    pivot_wider(names_from = sample_id, values_from = final_scale_intensity)
  
  # 移除包含"131"的列
  df_wide <- df_wide[, !grepl("131", names(df_wide))]
  
  # 合并Glycopeptides和Protein.Accession
  df_wide$Glycopeptides <- paste(df_wide$Glycopeptides, df_wide$Protein.Accession, sep = "-")
  df_wide <- df_wide %>% select(-Protein.Accession)
  
  # 替换列名为病人信息
  old_colnames <- colnames(df_wide)
  samples <- patient_info$Sample
  new_colnames <- patient_info$patients_repeat
  name_map <- setNames(new_colnames, samples)
  for (i in 1:length(old_colnames)) {
    if (old_colnames[i] %in% names(name_map)) {
      old_colnames[i] <- name_map[[old_colnames[i]]]
    }
  }
  colnames(df_wide) <- old_colnames
  
  # 检查列名是否完全对应
  unmatched <- setdiff(new_colnames, colnames(df_wide))
  if (length(unmatched) > 0) {
    stop(paste0("以下病人信息列未在数据中找到: ", paste(unmatched, collapse=", ")))
  }
  
  # 将NA替换为数值NA
  df_wide[df_wide == "NA"] <- NA
  
  # 计算有效值比例
  df_wide <- df_wide %>%
    rowwise() %>%
    mutate(valid_ratio = mean(!is.na(c_across(-Glycopeptides))))
  print("有效值比例统计：")
  print(summary(df_wide$valid_ratio))
  
  # 根据有效值比例筛选，若全丢失则自动降阈值
  df_filtered <- df_wide %>% filter(valid_ratio > valid_ratio) %>% select(-valid_ratio)
  if (nrow(df_filtered) == 0) {
    warning(sprintf("筛选后的糖肽数据为空，当前有效值比例阈值为%.2f，自动返回未筛选数据", valid_ratio))
    df_filtered <- df_wide %>% select(-valid_ratio)
  }
  
  # 使用KNN填补缺失值
  if (use_knn) {
    df_filled <- VIM::kNN(df_filtered, 
                        variable = setdiff(colnames(df_filtered), "Glycopeptides"),
                        k = knn_k)
    df_filled <- df_filled %>% select(-ends_with("_imp"))
  } else {
    df_filled <- df_filtered
  }
  
  return(df_filled)
}

# 计算样品重复相关性
calculate_sample_correlation <- function(data) {
  # 准备数据
  data_for_cor <- data
  rownames(data_for_cor) <- data_for_cor$Glycopeptides
  data_for_cor <- data_for_cor %>% select(-Glycopeptides)
  data_for_cor <- data_for_cor[, order(names(data_for_cor))]
  
  # 计算相关性矩阵
  correlation_matrix <- cor(data_for_cor, use = "pairwise.complete.obs")
  
  return(correlation_matrix)
}

# 计算样品均值
calculate_sample_means <- function(data) {
  # 转换为长格式
  df_long <- data %>%
    pivot_longer(
      cols = -Glycopeptides,
      names_to = c("Sample", "Repeat"),
      names_pattern = "(.*)_(repeat[0-9]+)"
    )
  
  # 计算每个样品的每个糖肽的均值
  df_mean <- df_long %>%
    group_by(Glycopeptides, Sample) %>%
    summarize(Mean_Expression = mean(value, na.rm = TRUE), .groups = 'drop')
  
  # 转换回宽格式
  df_wide_mean <- df_mean %>%
    pivot_wider(
      names_from = Sample,
      values_from = Mean_Expression
    )
  
  return(df_wide_mean)
}

# 执行Z-score标准化
perform_zscore_normalization <- function(data) {
  df_zscore <- data %>%
    mutate(across(where(is.numeric) & -Glycopeptides, ~ as.numeric(scale(.))))
  
  df_zscore$Glycopeptides <- data$Glycopeptides
  
  return(df_zscore)
}

# 执行PCA分析
perform_pca_analysis <- function(data) {
  # 提取糖肽名称
  gene_names <- data$Glycopeptides
  
  # 处理列名
  colnames(data)[-which(names(data) == "Glycopeptides")] <- 
    make.names(colnames(data)[-which(names(data) == "Glycopeptides")])
  
  # 提取注释信息
  annotations <- sapply(strsplit(colnames(data)[-which(names(data) == "Glycopeptides")], "\\."), 
                       function(x) x[length(x)])
  
  # 重命名列
  names(data)[-which(names(data) == "Glycopeptides")] <- 
    sapply(strsplit(colnames(data)[-which(names(data) == "Glycopeptides")], "\\."), 
           function(x) x[1])
  
  # 创建注释数据框
  annotation_df <- data.frame(
    Column = colnames(data)[-which(names(data) == "Glycopeptides")],
    Annotation = annotations
  )
  
  # 准备数值数据
  df_numeric <- data %>% select(-Glycopeptides)
  
  # 标准化
  df_scaled <- scale(df_numeric)
  
  # 执行PCA
  pca_result <- prcomp(df_scaled, center = TRUE, scale. = TRUE)
  
  # 准备结果
  pca_data <- as.data.frame(pca_result$x)
  pca_data$Sample <- rownames(pca_data)
  
  # 准备载荷数据
  loadings <- as.data.frame(pca_result$rotation)
  loadings$Annotation <- annotations
  
  return(list(
    pca_result = pca_result,
    pca_data = pca_data,
    loadings = loadings,
    annotation_df = annotation_df
  ))
}

# 绘制PCA图
plot_pca_results <- function(pca_results) {
  # 绘制载荷图
  p2 <- ggplot(pca_results$loadings, aes(x = PC1, y = PC2, color = Annotation)) +
    geom_point(size = 2) +
    labs(title = "PCA分析", 
         x = "主成分1载荷 (PC1)", 
         y = "主成分2载荷 (PC2)") +
    theme_classic()
  
  return(p2)
}

# 糖组数据可视化函数
plot_glycomics_data <- function(data_long, significant_results) {
  # 计算每个组的平均值
  data_summary <- data_long %>%
    group_by(Glyco_Molar_Mass, cancertype) %>%
    summarise(mean_intensity = mean(Value, na.rm = TRUE))
  
  # 按癌症类型归一化
  data_normalized <- data_summary %>%
    group_by(cancertype) %>%
    mutate(Mean_Value_normalized = (mean_intensity / max(mean_intensity)) * 100) %>%
    ungroup()
  
  # 获取显著性差异的糖链摩尔质量
  significant_glycos <- significant_results$Glyco_Molar_Mass
  
  # 绘制图表
  p <- ggplot(data_normalized, aes(x = Glyco_Molar_Mass, y = Mean_Value_normalized)) +
    geom_segment(aes(xend = Glyco_Molar_Mass, yend = 0), size = 1) +
    facet_wrap(~ cancertype, scales = "fixed", ncol = 1) +
    scale_x_continuous(limits = c(570, 5000)) +
    labs(x = "Glyco_Molar_Mass", y = "Mean_Value") +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    geom_text(data = subset(data_normalized, Glyco_Molar_Mass %in% significant_glycos),
              aes(label = Glyco_Molar_Mass), color = "red", size = 3, vjust = -0.5)
  
  return(p)
}

# 自动保存CSV函数
auto_save_csv <- function(df, filename) {
  write.csv(df, filename, row.names = FALSE)
}

# IGP定量分析与热图绘制函数
igp_quant_heatmap_analysis <- function(expr_df, group_df, ratio_upper = 1.2, ratio_lower = 0.83, color_palette = colorRampPalette(c("blue", "white", "red"))(50), breaks = seq(-4, 4, length.out = 50)) {
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  # 动态提取所有实际存在的分组类型
  group_levels <- unique(group_df$group)
  colnames(expr_df) <- as.character(expr_df[1,])
  expr_df <- expr_df[-1,]
  gene_expression_long <- tidyr::gather(expr_df, key = "sample", value = "expression", -Glycopeptides)
  merged_data <- merge(gene_expression_long, group_df, by.x = "sample", by.y = "sample")
  wide_data <- tidyr::spread(merged_data, key = Glycopeptides, value = expression)
  wide_data[ , 3:ncol(wide_data)] <- lapply(wide_data[ , 3:ncol(wide_data)], as.numeric)
  data_long <- wide_data %>%
    tidyr::gather(key = "gene", value = "expression", -sample, -group)
  data_long$expression <- as.numeric(data_long$expression)
  # 后续所有分组相关分析都用group_levels
  # ... existing code ...
  # 返回时加上group_levels
  return(list(
    crc_heatmap = crc_heatmap,
    lihc_heatmap = lihc_heatmap,
    combined_heatmap = combined_heatmap,
    combined_heatmap_nocluster = combined_heatmap_nocluster,
    zscore_heatmap = zscore_heatmap,
    zscore_heatmap_cluster = zscore_heatmap_cluster,
    zscore_heatmap_nocluster = zscore_heatmap_nocluster,
    group_levels = group_levels
  ))
}

# 区域热图分析
regional_analysis_and_plots <- function(expr_df, group_df, glycan_type_df, color_palette = colorRampPalette(c("blue", "white", "red"))(50), breaks = seq(-4, 4, length.out = 50)) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  
  # --- 数据预处理 ---
  colnames(expr_df) <- as.character(expr_df[1,])
  expr_df <- expr_df[-1,]
  gene_expression_long <- tidyr::gather(expr_df, key = "sample", value = "expression", -Glycopeptides)
  merged_data <- merge(gene_expression_long, group_df, by.x = "sample", by.y = "sample")
  wide_data <- tidyr::spread(merged_data, key = Glycopeptides, value = expression)
  wide_data[ , 3:ncol(wide_data)] <- lapply(wide_data[ , 3:ncol(wide_data)], as.numeric)
  
  # --- 计算均值比 ---
  data_long <- wide_data %>%
    gather(key = "gene", value = "expression", -sample, -group)
  data_long$expression <- as.numeric(data_long$expression)
  group_means <- data_long %>%
    group_by(group, gene) %>%
    summarize(mean_expression = mean(expression, na.rm = TRUE), .groups = 'drop')
  group_means_wide <- group_means %>%
    spread(key = group, value = mean_expression)
    
  differences <- group_means_wide %>%
    mutate(CRLM_minus_CRC = CRLM / CRC,
           CRLM_minus_LIHC = CRLM / LIHC) %>%
    dplyr::select(gene, CRLM_minus_CRC, CRLM_minus_LIHC)
    
  # --- 生成散点图 ---
  differences_for_plot <- differences
  differences_for_plot <- tidyr::separate(differences_for_plot, gene, sep = "-", into=c("peptides","glycan","Protein"), remove = FALSE)
  differences_for_plot <- merge(differences_for_plot, glycan_type_df)
  glycan_colors <- c("High-mannose" = "#AFC778",
                     "Fucosylation" = "#9C9BE9",
                     "Paucimannose" = "#A3CDEA",
                     "Sialytation" = "#F09496",
                     "Fuc-Sia" = "#F4DFDD",
                     "Others"="#EFDBB9")
  
  scatter_plot <- ggplot(differences_for_plot, aes(x = CRLM_minus_LIHC, y = CRLM_minus_CRC, color = glycan.type)) +
    geom_point() +
    geom_rect(aes(xmin = 0.83, xmax = 1.2, ymin = 0.83, ymax = 1.2), color = "red", fill = NA) +
    geom_hline(yintercept = 0.83, linetype = "solid", color = "green") +
    geom_hline(yintercept = 1.2, linetype = "solid", color = "green") +
    geom_vline(xintercept = 0.83, linetype = "solid", color = "green") +
    geom_vline(xintercept = 1.2, linetype = "solid", color = "green") +
    scale_color_manual(values = glycan_colors) +
    labs(title = "Scatter plot of gene expression differences",
         x = "CRLM / LIHC",
         y = "CRLM / CRC") +
    theme_classic()
    
  # --- 区域划分 ---
  regions_df <- differences %>%
    mutate(region = case_when(
      CRLM_minus_LIHC < 0.83 & CRLM_minus_CRC < 0.83 ~ "Bottom-Left",
      CRLM_minus_LIHC < 0.83 & CRLM_minus_CRC >= 0.83 & CRLM_minus_CRC <= 1.2 ~ "Left-Center",
      CRLM_minus_LIHC < 0.83 & CRLM_minus_CRC > 1.2 ~ "Top-Left",
      CRLM_minus_LIHC >= 0.83 & CRLM_minus_LIHC <= 1.2 & CRLM_minus_CRC < 0.83 ~ "Bottom-Center",
      CRLM_minus_LIHC >= 0.83 & CRLM_minus_LIHC <= 1.2 & CRLM_minus_CRC >= 0.83 & CRLM_minus_CRC <= 1.2 ~ "Center",
      CRLM_minus_LIHC >= 0.83 & CRLM_minus_LIHC <= 1.2 & CRLM_minus_CRC > 1.2 ~ "Top-Center",
      CRLM_minus_LIHC > 1.2 & CRLM_minus_CRC < 0.83 ~ "Bottom-Right",
      CRLM_minus_LIHC > 1.2 & CRLM_minus_CRC >= 0.83 & CRLM_minus_CRC <= 1.2 ~ "Right-Center",
      CRLM_minus_LIHC > 1.2 & CRLM_minus_CRC > 1.2 ~ "Top-Right"
    ))
  region_counts <- as.data.frame(table(regions_df$region))
  colnames(region_counts) <- c("Region", "Count")
  
  # --- 生成热图 ---
  # 准备差异数据
  crlm_means <- data_long %>%
    filter(group == "CRLM") %>%
    group_by(gene) %>%
    summarize(mean_expression = mean(expression, na.rm = TRUE))
  crc_data_long <- data_long %>% filter(group == "CRC")
  lihc_data_long <- data_long %>% filter(group == "LIHC")
  crc_vs_crlm <- crc_data_long %>%
    left_join(crlm_means, by = "gene") %>%
    mutate(expression_diff = mean_expression / expression) %>%
    dplyr::select(sample, gene, expression_diff)
  crc_vs_crlm_wide <- crc_vs_crlm %>% spread(key = gene, value = expression_diff)
  lihc_vs_crlm <- lihc_data_long %>%
    left_join(crlm_means, by = "gene") %>%
    mutate(expression_diff = mean_expression / expression) %>%
    dplyr::select(sample, gene, expression_diff)
  lihc_vs_crlm_wide <- lihc_vs_crlm %>% spread(key = gene, value = expression_diff)
  
  # 合并差异数据
  combined_wide_diff <- rbind(crc_vs_crlm_wide, lihc_vs_crlm_wide)
  genes_not_center <- regions_df %>% filter(region != "Center") %>% pull(gene)
  
  # 差异热图
  heatmap_diff_data <- combined_wide_diff[, c("sample", genes_not_center)]
  rownames(heatmap_diff_data) <- heatmap_diff_data$sample
  heatmap_diff_data <- heatmap_diff_data[, -1] 
  heatmap_diff_data <- as.matrix(heatmap_diff_data)
  
  heatmap_diff <- pheatmap(
    t(heatmap_diff_data),
    cluster_rows = TRUE, 
    cluster_cols = FALSE, 
    show_colnames = TRUE,
    show_rownames = TRUE,
    color = color_palette,
    breaks = breaks,
    main = "Heatmap of genes not in Center region (Expression Difference)"
  )
  
  # z-score热图
  zscore_heatmap_data <- wide_data %>%
    dplyr::select(sample, all_of(genes_not_center))
  col_annotation_df <- wide_data[, c("sample", "group")]
  rownames(col_annotation_df) <- col_annotation_df$sample
  col_annotation_df <- col_annotation_df %>% dplyr::select(group)
  colnames(col_annotation_df) <- "Cancer type"
  
  rownames(zscore_heatmap_data) <- zscore_heatmap_data$sample
  zscore_heatmap_data <- zscore_heatmap_data[, -1] 
  
  ordered_samples <- rownames(col_annotation_df)[order(col_annotation_df$`Cancer type`)]
  zscore_heatmap_data <- zscore_heatmap_data[ordered_samples, ]
  col_annotation_df <- col_annotation_df[ordered_samples, , drop = FALSE]
  
  heatmap_zscore <- pheatmap(
    t(zscore_heatmap_data),
    cluster_rows = TRUE, 
    cluster_cols = FALSE, 
    show_colnames = TRUE,
    show_rownames = TRUE,
    scale = "row", 
    color = color_palette,
    breaks = breaks,
    annotation_col = col_annotation_df,
    main = "Heatmap of genes not in Center region (Z-score)"
  )
  
  return(list(
    scatter_plot = scatter_plot,
    region_counts = region_counts,
    heatmap_diff = heatmap_diff,
    heatmap_zscore = heatmap_zscore
  ))
}

# 韦恩图分析函数
venn_diagram_analysis <- function(expr_df, group_df, glycan_type_df) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggvenn)
  library(ggVennDiagram)
  library(patchwork)
  
  # --- 数据预处理 ---
  colnames(expr_df) <- as.character(expr_df[1,])
  expr_df <- expr_df[-1,]
  gene_expression_long <- tidyr::gather(expr_df, key = "sample", value = "expression", -Glycopeptides)
  merged_data <- merge(gene_expression_long, group_df, by.x = "sample", by.y = "sample")
  wide_data <- tidyr::spread(merged_data, key = Glycopeptides, value = expression)
  wide_data[ , 3:ncol(wide_data)] <- lapply(wide_data[ , 3:ncol(wide_data)], as.numeric)
  
  # --- 计算均值比 ---
  data_long <- wide_data %>%
    gather(key = "gene", value = "expression", -sample, -group)
  data_long$expression <- as.numeric(data_long$expression)
  group_means <- data_long %>%
    group_by(group, gene) %>%
    summarize(mean_expression = mean(expression, na.rm = TRUE), .groups = 'drop')
  group_means_wide <- group_means %>%
    spread(key = group, value = mean_expression)
    
  differences <- group_means_wide %>%
    mutate(CRLM_minus_CRC = CRLM / CRC,
           CRLM_minus_LIHC = CRLM / LIHC) %>%
    dplyr::select(gene, CRLM_minus_CRC, CRLM_minus_LIHC)
    
  # --- 区域划分 ---
  regions_df <- differences %>%
    mutate(region = case_when(
      CRLM_minus_LIHC < 0.83 & CRLM_minus_CRC < 0.83 ~ "Bottom-Left",
      CRLM_minus_LIHC < 0.83 & CRLM_minus_CRC >= 0.83 & CRLM_minus_CRC <= 1.2 ~ "Left-Center",
      CRLM_minus_LIHC < 0.83 & CRLM_minus_CRC > 1.2 ~ "Top-Left",
      CRLM_minus_LIHC >= 0.83 & CRLM_minus_LIHC <= 1.2 & CRLM_minus_CRC < 0.83 ~ "Bottom-Center",
      CRLM_minus_LIHC >= 0.83 & CRLM_minus_LIHC <= 1.2 & CRLM_minus_CRC >= 0.83 & CRLM_minus_CRC <= 1.2 ~ "Center",
      CRLM_minus_LIHC >= 0.83 & CRLM_minus_LIHC <= 1.2 & CRLM_minus_CRC > 1.2 ~ "Top-Center",
      CRLM_minus_LIHC > 1.2 & CRLM_minus_CRC < 0.83 ~ "Bottom-Right",
      CRLM_minus_LIHC > 1.2 & CRLM_minus_CRC >= 0.83 & CRLM_minus_CRC <= 1.2 ~ "Right-Center",
      CRLM_minus_LIHC > 1.2 & CRLM_minus_CRC > 1.2 ~ "Top-Right"
    ))
  
  # 获取非中心区域的基因
  different_gene_mean_calculate_regions <- regions_df %>% 
    filter(region != "Center") %>% 
    pull(gene)
  
  # --- 热图分析获取combined_genes ---
  # 准备差异数据
  crlm_means <- data_long %>%
    filter(group == "CRLM") %>%
    group_by(gene) %>%
    summarize(mean_expression = mean(expression, na.rm = TRUE))
  
  # CRC vs CRLM
  crc_data <- wide_data %>% filter(group == "CRC")
  crc_data_long <- crc_data %>% gather(key = "gene", value = "expression", -sample, -group)
  crc_data_long$expression <- as.numeric(crc_data_long$expression)
  crc_vs_crlm <- crc_data_long %>%
    left_join(crlm_means, by = "gene") %>%
    mutate(expression_diff = mean_expression / expression) %>%
    dplyr::select(sample, gene, expression_diff)
  crc_vs_crlm_wide <- crc_vs_crlm %>% spread(key = gene, value = expression_diff)
  crc_expression_diff_matrix <- as.matrix(crc_vs_crlm_wide[,-1])
  rownames(crc_expression_diff_matrix) <- crc_vs_crlm_wide$sample
  crc_gene_filter <- apply(crc_expression_diff_matrix, 2, function(x) {
    sum(x > 1.2) >= length(x)/2 || sum(x < 0.83) >= length(x)/2
  })
  crc_filtered_genes <- colnames(crc_expression_diff_matrix)[crc_gene_filter]
  
  # LIHC vs CRLM
  lihc_data <- wide_data %>% filter(group == "LIHC")
  lihc_data_long <- lihc_data %>% gather(key = "gene", value = "expression", -sample, -group)
  lihc_data_long$expression <- as.numeric(lihc_data_long$expression)
  lihc_vs_crlm <- lihc_data_long %>%
    left_join(crlm_means, by = "gene") %>%
    mutate(expression_diff = mean_expression / expression) %>%
    dplyr::select(sample, gene, expression_diff)
  lihc_vs_crlm_wide <- lihc_vs_crlm %>% spread(key = gene, value = expression_diff)
  lihc_expression_diff_matrix <- as.matrix(lihc_vs_crlm_wide[,-1])
  rownames(lihc_expression_diff_matrix) <- lihc_vs_crlm_wide$sample
  lihc_gene_filter <- apply(lihc_expression_diff_matrix, 2, function(x) {
    sum(x > 1.2) >= length(x)/2 || sum(x < 0.83) >= length(x)/2
  })
  lihc_filtered_genes <- colnames(lihc_expression_diff_matrix)[lihc_gene_filter]
  
  # 合并筛选基因
  combined_genes <- union(crc_filtered_genes, lihc_filtered_genes)
  
  # --- 韦恩图1：两种方法确定的糖肽 ---
  venn_data_1 <- list(
    Region = different_gene_mean_calculate_regions,
    Combined = combined_genes
  )
  
  venn_plot_1 <- ggvenn(venn_data_1, 
                        fill_color = c("blue", "red"),
                        stroke_size = 0.5, 
                        set_name_size = 4,
                        text_size = 5) +
    labs(title = "两种方法确定的糖肽") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  # --- 统计检验获取显著性差异基因 ---
  different_genes_from_both <- union(different_gene_mean_calculate_regions, combined_genes)
  
  # 初始化存储显著性差异基因的向量
  significant_genes <- c()
  
  # 遍历每个基因，进行t检验
  for (gene in different_genes_from_both) {
    gene_expression <- wide_data[, c("sample", "group", gene)]
    colnames(gene_expression) <- c("sample", "group", "expression")
    
    # 对CRC和CRLM进行t检验
    t_test_crc_crlm <- t.test(expression ~ group, data = subset(gene_expression, group %in% c("CRC", "CRLM")))
    
    # 对CRLM和LIHC进行t检验
    t_test_crlm_lihc <- t.test(expression ~ group, data = subset(gene_expression, group %in% c("CRLM", "LIHC")))
    
    # 检查p值是否小于0.05
    if (t_test_crc_crlm$p.value < 0.05 || t_test_crlm_lihc$p.value < 0.05) {
      significant_genes <- c(significant_genes, gene)
    }
  }
  
  significant_genes <- unique(significant_genes)
  
  # 分别筛选CRC和LIHC的显著性差异基因
  final_significant_genes_crc <- c()
  final_significant_genes_lihc <- c()
  
  for (gene in significant_genes) {
    gene_expression <- wide_data[, c("sample", "group", gene)]
    colnames(gene_expression) <- c("sample", "group", "expression")
    
    t_test_crc_crlm <- t.test(expression ~ group, data = subset(gene_expression, group %in% c("CRC", "CRLM")))
    t_test_crlm_lihc <- t.test(expression ~ group, data = subset(gene_expression, group %in% c("CRLM", "LIHC")))
    
    if (t_test_crc_crlm$p.value < 0.05) {
      final_significant_genes_crc <- c(final_significant_genes_crc, gene)
    }
    if (t_test_crlm_lihc$p.value < 0.05) {
      final_significant_genes_lihc <- c(final_significant_genes_lihc, gene)
    }
  }
  
  final_significant_genes_crc <- unique(final_significant_genes_crc)
  final_significant_genes_lihc <- unique(final_significant_genes_lihc)
  
  # --- 韦恩图2：显著性差异基因 ---
  venn_data_2 <- list(
    All_sig = significant_genes,
    CRLM_vs_CRC = final_significant_genes_crc,
    CRLM_vs_LIHC = final_significant_genes_lihc
  )
  
  venn_plot_2 <- ggvenn(venn_data_2, 
                        fill_color = c("blue", "red", "pink"),
                        stroke_size = 0.5, 
                        set_name_size = 4,
                        text_size = 5) +
    labs(title = "Venn Diagram of Significant Genes") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = "none"
    )
  
  # --- 韦恩图3：使用ggVennDiagram ---
  venn_plot_3 <- ggVennDiagram(venn_data_2, 
                              fill = c("blue", "red", "pink"), 
                              label_alpha = 0, 
                              label_size = 4) +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(title = "Venn Diagram of Significant Genes") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.position = "none"
    )
  
  # --- 显著性差异基因热图 ---
  wide_data_for_hp <- wide_data %>%
    dplyr::select(significant_genes)
  rownames(wide_data_for_hp) <- wide_data$sample
  
  # 行注释
  row_annotation <- as.data.frame(colnames(wide_data_for_hp))
  rownames(row_annotation) <- row_annotation$`colnames(wide_data_for_hp)`
  colnames(row_annotation) <- "Glycopeptides"
  row_annotation$Glycopeptides1 <- row_annotation$Glycopeptides
  row_annotation <- separate(row_annotation, Glycopeptides, sep = "-", into = c("Peptides", "glycan", "Proteins"))
  row_annotation <- merge(row_annotation, glycan_type_df)
  a <- row_annotation$Glycopeptides1
  rownames(row_annotation) <- a
  row_annotation <- row_annotation[, -c(1, 2, 4)]
  
  wide_data_for_hp <- t(wide_data_for_hp)
  
  # 列注释
  col_annotation <- wide_data[, c(1, 2)]
  col_annotation <- as.data.frame(col_annotation[, -1], drop = FALSE)
  rownames(col_annotation) <- wide_data$sample
  colnames(col_annotation) <- "Cancer type"
  wide_data_for_hp <- wide_data_for_hp[, order(col_annotation$`Cancer type`)]
  
  # 定义糖链类型颜色
  glycan_colors <- c("High-mannose" = "#AFC778",
                     "Fucosylation" = "#9C9BE9",
                     "Paucimannose" = "#A3CDEA",
                     "Sialytation" = "#F09496",
                     "Fuc-Sia" = "#F4DFDD",
                     "Others" = "#EFDBB9")
  
  # 显著性差异基因热图
  significant_heatmap <- pheatmap(wide_data_for_hp,
                                 cluster_rows = TRUE, 
                                 cluster_cols = FALSE, 
                                 show_colnames = TRUE,
                                 show_rownames = TRUE,
                                 scale = "row", 
                                 color = colorRampPalette(c("blue", "white", "red"))(50),
                                 breaks = seq(-4, 4, length.out = 50),
                                 annotation_col = col_annotation,
                                 annotation_row = row_annotation,
                                 main = "Heatmap of 所有筛选出来具有显著性差异的",
                                 annotation_colors = list(glycan.type = glycan_colors))
  
  # --- 仅CRLM样本的显著性差异基因热图 ---
  wide_data_only_CRLM <- subset(wide_data, wide_data$group == "CRLM")
  wide_data_for_hp_crlm <- wide_data_only_CRLM %>%
    dplyr::select(significant_genes)
  rownames(wide_data_for_hp_crlm) <- wide_data_only_CRLM$sample
  
  # 行注释（与上面相同）
  row_annotation_crlm <- as.data.frame(colnames(wide_data_for_hp_crlm))
  rownames(row_annotation_crlm) <- row_annotation_crlm$`colnames(wide_data_for_hp_crlm)`
  colnames(row_annotation_crlm) <- "Glycopeptides"
  row_annotation_crlm$Glycopeptides1 <- row_annotation_crlm$Glycopeptides
  row_annotation_crlm <- separate(row_annotation_crlm, Glycopeptides, sep = "-", into = c("Peptides", "glycan", "Proteins"))
  row_annotation_crlm <- merge(row_annotation_crlm, glycan_type_df)
  a_crlm <- row_annotation_crlm$Glycopeptides1
  rownames(row_annotation_crlm) <- a_crlm
  row_annotation_crlm <- row_annotation_crlm[, -c(1, 2, 4)]
  
  wide_data_for_hp_crlm <- t(wide_data_for_hp_crlm)
  
  # 列注释
  col_annotation_crlm <- wide_data_only_CRLM[, c(1, 2)]
  col_annotation_crlm <- as.data.frame(col_annotation_crlm[, -1], drop = FALSE)
  rownames(col_annotation_crlm) <- wide_data_only_CRLM$sample
  colnames(col_annotation_crlm) <- "Cancer type"
  wide_data_for_hp_crlm <- wide_data_for_hp_crlm[, order(col_annotation_crlm$`Cancer type`)]
  
  # CRLM显著性差异基因热图
  crlm_significant_heatmap <- pheatmap(wide_data_for_hp_crlm,
                                      cluster_rows = TRUE, 
                                      cluster_cols = TRUE, 
                                      show_colnames = TRUE,
                                      show_rownames = TRUE,
                                      scale = "row", 
                                      color = colorRampPalette(c("blue", "white", "red"))(50),
                                      breaks = seq(-4, 4, length.out = 50),
                                      annotation_col = col_annotation_crlm,
                                      annotation_row = row_annotation_crlm,
                                      main = "Heatmap of 所有筛选出来具有显著性差异的 (CRLM only)",
                                      annotation_colors = list(glycan.type = glycan_colors))
  
  # --- 按糖链类型排序的CRLM热图 ---
  row_annotation_crlm$Glycan_Type_Factor <- factor(row_annotation_crlm$glycan.type, 
                                                  levels = c("High-mannose", "Fucosylation", "Paucimannose", "Sialytation", "Fuc-Sia", "Others"))
  row_annotation_crlm <- row_annotation_crlm[order(row_annotation_crlm$Glycan_Type_Factor), ]
  wide_data_for_hp_sorted <- wide_data_for_hp_crlm[rownames(row_annotation_crlm), ]
  
  crlm_sorted_heatmap <- pheatmap(wide_data_for_hp_sorted,
                                 cluster_rows = FALSE,  # 关闭行聚类，因为我们已经手动排序
                                 cluster_cols = TRUE, 
                                 show_colnames = TRUE,
                                 show_rownames = TRUE,
                                 scale = "row", 
                                 color = colorRampPalette(c("blue", "white", "red"))(50),
                                 breaks = seq(-4, 4, length.out = 50),
                                 annotation_col = col_annotation_crlm,
                                 annotation_row = row_annotation_crlm,
                                 main = "Heatmap of 所有筛选出来具有显著性差异的 (按糖链类型排序)",
                                 annotation_colors = list(glycan.type = glycan_colors))
  
  # 组合3个venn图
  venn_combined_plot <- (venn_plot_1 + venn_plot_2 + venn_plot_3) +
    plot_layout(ncol = 3, guides = "collect") +
    plot_annotation(
      title = "Venn Diagram Analysis Results",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
  
  return(list(
    venn_plot_1 = venn_plot_1,
    venn_plot_2 = venn_plot_2,
    venn_plot_3 = venn_plot_3,
    significant_heatmap = significant_heatmap,
    crlm_significant_heatmap = crlm_significant_heatmap,
    crlm_sorted_heatmap = crlm_sorted_heatmap,
    significant_genes = significant_genes,
    final_significant_genes_crc = final_significant_genes_crc,
    final_significant_genes_lihc = final_significant_genes_lihc,
    different_gene_mean_calculate_regions = different_gene_mean_calculate_regions,
    combined_genes = combined_genes,
    venn_combined_plot = venn_combined_plot
  ))
}

# LDA图和GO分析函数
lda_go_analysis <- function(expr_df, group_df, glycan_type_df) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(factoextra)
  library(MASS)
  library(e1071)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  # --- 数据预处理 ---
  colnames(expr_df) <- as.character(expr_df[1,])
  expr_df <- expr_df[-1,]
  gene_expression_long <- tidyr::gather(expr_df, key = "sample", value = "expression", -Glycopeptides)
  merged_data <- merge(gene_expression_long, group_df, by.x = "sample", by.y = "sample")
  wide_data <- tidyr::spread(merged_data, key = Glycopeptides, value = expression)
  wide_data[ , 3:ncol(wide_data)] <- lapply(wide_data[ , 3:ncol(wide_data)], as.numeric)
  
  # --- 计算均值比 ---
  data_long <- wide_data %>%
    gather(key = "gene", value = "expression", -sample, -group)
  data_long$expression <- as.numeric(data_long$expression)
  group_means <- data_long %>%
    group_by(group, gene) %>%
    summarize(mean_expression = mean(expression, na.rm = TRUE), .groups = 'drop')
  group_means_wide <- group_means %>%
    spread(key = group, value = mean_expression)
    
  differences <- group_means_wide %>%
    mutate(CRLM_minus_CRC = CRLM / CRC,
           CRLM_minus_LIHC = CRLM / LIHC) %>%
    dplyr::select(gene, CRLM_minus_CRC, CRLM_minus_LIHC)
    
  # --- 区域划分 ---
  regions_df <- differences %>%
    mutate(region = case_when(
      CRLM_minus_LIHC < 0.83 & CRLM_minus_CRC < 0.83 ~ "Bottom-Left",
      CRLM_minus_LIHC < 0.83 & CRLM_minus_CRC >= 0.83 & CRLM_minus_CRC <= 1.2 ~ "Left-Center",
      CRLM_minus_LIHC < 0.83 & CRLM_minus_CRC > 1.2 ~ "Top-Left",
      CRLM_minus_LIHC >= 0.83 & CRLM_minus_LIHC <= 1.2 & CRLM_minus_CRC < 0.83 ~ "Bottom-Center",
      CRLM_minus_LIHC >= 0.83 & CRLM_minus_LIHC <= 1.2 & CRLM_minus_CRC >= 0.83 & CRLM_minus_CRC <= 1.2 ~ "Center",
      CRLM_minus_LIHC >= 0.83 & CRLM_minus_LIHC <= 1.2 & CRLM_minus_CRC > 1.2 ~ "Top-Center",
      CRLM_minus_LIHC > 1.2 & CRLM_minus_CRC < 0.83 ~ "Bottom-Right",
      CRLM_minus_LIHC > 1.2 & CRLM_minus_CRC >= 0.83 & CRLM_minus_CRC <= 1.2 ~ "Right-Center",
      CRLM_minus_LIHC > 1.2 & CRLM_minus_CRC > 1.2 ~ "Top-Right"
    ))
  
  # 获取非中心区域的基因
  different_gene_mean_calculate_regions <- regions_df %>% 
    filter(region != "Center") %>% 
    pull(gene)
  
  # --- 热图分析获取combined_genes ---
  # 准备差异数据
  crlm_means <- data_long %>%
    filter(group == "CRLM") %>%
    group_by(gene) %>%
    summarize(mean_expression = mean(expression, na.rm = TRUE))
  
  # CRC vs CRLM
  crc_data <- wide_data %>% filter(group == "CRC")
  crc_data_long <- crc_data %>% gather(key = "gene", value = "expression", -sample, -group)
  crc_data_long$expression <- as.numeric(crc_data_long$expression)
  crc_vs_crlm <- crc_data_long %>%
    left_join(crlm_means, by = "gene") %>%
    mutate(expression_diff = mean_expression / expression) %>%
    dplyr::select(sample, gene, expression_diff)
  crc_vs_crlm_wide <- crc_vs_crlm %>% spread(key = gene, value = expression_diff)
  crc_expression_diff_matrix <- as.matrix(crc_vs_crlm_wide[,-1])
  rownames(crc_expression_diff_matrix) <- crc_vs_crlm_wide$sample
  crc_gene_filter <- apply(crc_expression_diff_matrix, 2, function(x) {
    sum(x > 1.2) >= length(x)/2 || sum(x < 0.83) >= length(x)/2
  })
  crc_filtered_genes <- colnames(crc_expression_diff_matrix)[crc_gene_filter]
  
  # LIHC vs CRLM
  lihc_data <- wide_data %>% filter(group == "LIHC")
  lihc_data_long <- lihc_data %>% gather(key = "gene", value = "expression", -sample, -group)
  lihc_data_long$expression <- as.numeric(lihc_data_long$expression)
  lihc_vs_crlm <- lihc_data_long %>%
    left_join(crlm_means, by = "gene") %>%
    mutate(expression_diff = mean_expression / expression) %>%
    dplyr::select(sample, gene, expression_diff)
  lihc_vs_crlm_wide <- lihc_vs_crlm %>% spread(key = gene, value = expression_diff)
  lihc_expression_diff_matrix <- as.matrix(lihc_vs_crlm_wide[,-1])
  rownames(lihc_expression_diff_matrix) <- lihc_vs_crlm_wide$sample
  lihc_gene_filter <- apply(lihc_expression_diff_matrix, 2, function(x) {
    sum(x > 1.2) >= length(x)/2 || sum(x < 0.83) >= length(x)/2
  })
  lihc_filtered_genes <- colnames(lihc_expression_diff_matrix)[lihc_gene_filter]
  
  # 合并筛选基因
  combined_genes <- union(crc_filtered_genes, lihc_filtered_genes)
  
  # --- 统计检验获取显著性差异基因 ---
  different_genes_from_both <- union(different_gene_mean_calculate_regions, combined_genes)
  
  # 初始化存储显著性差异基因的向量
  significant_genes <- c()
  
  # 遍历每个基因，进行t检验
  for (gene in different_genes_from_both) {
    gene_expression <- wide_data[, c("sample", "group", gene)]
    colnames(gene_expression) <- c("sample", "group", "expression")
    
    # 对CRC和CRLM进行t检验
    t_test_crc_crlm <- t.test(expression ~ group, data = subset(gene_expression, group %in% c("CRC", "CRLM")))
    
    # 对CRLM和LIHC进行t检验
    t_test_crlm_lihc <- t.test(expression ~ group, data = subset(gene_expression, group %in% c("CRLM", "LIHC")))
    
    # 检查p值是否小于0.05
    if (t_test_crc_crlm$p.value < 0.05 || t_test_crlm_lihc$p.value < 0.05) {
      significant_genes <- c(significant_genes, gene)
    }
  }
  
  significant_genes <- unique(significant_genes)
  
  # --- LDA分析 ---
  # 准备基因表达数据
  gene_expression <- wide_data[, c("sample", "group", significant_genes)]
  gene_expression$sample <- as.factor(gene_expression$sample)
  gene_expression$group <- as.factor(gene_expression$group)
  
  # 提取基因表达数据部分
  gene_data <- gene_expression[, 3:ncol(gene_expression)]
  
  # 进行PCA分析
  pca_result <- prcomp(gene_data, scale. = TRUE)
  
  # 创建PCA图
  pca_plot <- fviz_pca_ind(pca_result, 
                           geom.ind = "point", 
                           col.ind = gene_expression$group, 
                           addEllipses = TRUE, 
                           legend.title = "Group") +
    labs(title = "PCA Analysis of Significant Genes") +
    theme_classic()
  
  # 进行LDA分析
  lda_result <- lda(group ~ ., data = gene_expression[, -1])
  
  # 预测并可视化LDA结果
  lda_values <- predict(lda_result)
  lda_df <- data.frame(lda_values$x, group = gene_expression$group)
  
  # 可视化LDA结果
  lda_plot <- ggplot(lda_df, aes(x = LD1, y = LD2, color = group)) +
    geom_point() +
    stat_ellipse(aes(fill = group), alpha = 0.2) +
    labs(title = "LDA - Linear Discriminant Analysis") +
    theme_classic()
  
  # --- SVM分析 ---
  # 提取前两个主成分
  pca_data <- as.data.frame(pca_result$x[, 1:2])
  colnames(pca_data) <- c("PC1", "PC2")
  
  # 将PCA结果与原数据框合并
  gene_expression_pca <- cbind(gene_expression, pca_data)
  
  # 进行SVM分析
  svm_model <- svm(group ~ ., data = gene_expression_pca[, c("group", "PC1", "PC2")], kernel = "linear")
  
  # 预测结果
  svm_prediction <- predict(svm_model, gene_expression_pca[, c("PC1", "PC2")])
  gene_expression_pca$svm_prediction <- svm_prediction
  
  # 可视化SVM结果
  svm_plot <- ggplot(gene_expression_pca, aes(x = PC1, y = PC2, color = svm_prediction)) +
    geom_point() +
    labs(title = "SVM - Support Vector Machine") +
    theme_classic()
  
  # --- GO分析 ---
  # 准备显著性差异基因数据
  significant_genes_for_GO_analysis <- as.data.frame(significant_genes)
  significant_genes_for_GO_analysis <- significant_genes_for_GO_analysis %>%
    separate(significant_genes, into = c("Peptides", "glycan", "Proteins"), sep = "-")
  unique_proteins <- unique(significant_genes_for_GO_analysis$Proteins)
  # 将UniProt ID转换为Entrez ID
  gene_ids <- tryCatch({
    bitr(unique_proteins, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) NULL)
  # 进行GO富集分析
  go_results <- tryCatch({
    if (!is.null(gene_ids) && nrow(gene_ids) > 0) {
      enrichGO(gene = gene_ids$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "BP",
               pAdjustMethod = "BH",
               qvalueCutoff = 0.05)
    } else {
      NULL
    }
  }, error = function(e) NULL)
  if (is.null(go_results) || nrow(summary(go_results)) == 0) {
    go_barplot <- NULL
    go_dotplot <- NULL
    enrichment_results <- NULL
  } else {
    go_barplot <- barplot(go_results, showCategory = 20) +
      labs(title = "GO Enrichment Analysis - Bar Plot") +
      theme_classic()
    go_dotplot <- dotplot(go_results, showCategory = 20) +
      labs(title = "GO Enrichment Analysis - Dot Plot") +
      theme_classic()
    # ... enrichment_results 相关代码 ...
    go_summary <- summary(go_results)
    go_terms <- go_summary[, c("ID", "Description")]
    go_proteins <- go_results@result[, c("ID", "geneID")]
    enrichment_results <- merge(go_terms, go_proteins, by.x = "ID", by.y = "ID")
    enrichment_results$Proteins <- sapply(enrichment_results$geneID, function(x) paste(unlist(strsplit(as.character(x), "/")), collapse = ", "))
    enrichment_results <- enrichment_results[, c("ID", "Description", "Proteins")]
    convert_to_gene_symbols <- function(entrez_ids) {
      gene_symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
      return(gene_symbols)
    }
    enrichment_results$Proteins_GeneSymbols <- sapply(strsplit(as.character(enrichment_results$Proteins), ", "), function(ids) {
      gene_symbols <- convert_to_gene_symbols(ids)
      paste(gene_symbols, collapse = ", ")
    })
  }
  # ... existing code ...
  
  # --- 糖链类型分析 ---
  # 合并糖链类型信息
  significant_genes_for_GO_analysis <- merge(significant_genes_for_GO_analysis, glycan_type_df)
  
  # 统计糖链类型分布
  glycan_type_counts <- table(significant_genes_for_GO_analysis$glycan.type)
  glycan_data <- data.frame(
    glycan_type = names(glycan_type_counts),
    count = as.numeric(glycan_type_counts)
  )
  
  # 糖链类型饼图
  glycan_pie_plot <- ggplot(glycan_data, aes(x = "", y = count, fill = glycan_type)) +
    geom_bar(stat = "identity", width = 1) + 
    coord_polar(theta = "y") +
    labs(title = "Glycan Type Distribution") +
    theme_void() +
    theme(legend.title = element_blank()) +
    geom_text(aes(label = paste0(glycan_type, "\n", round(count/sum(count)*100, 1), "%")),
              position = position_stack(vjust = 0.5))
  
  # 糖链类型柱状图
  glycan_bar_plot <- ggplot(glycan_data, aes(x = glycan_type, y = count, fill = glycan_type)) +
    geom_bar(stat = "identity") +
    labs(title = "Glycan Type Distribution", x = "Glycan Type", y = "Count") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(
    pca_plot = pca_plot,
    lda_plot = lda_plot,
    svm_plot = svm_plot,
    go_barplot = go_barplot,
    go_dotplot = go_dotplot,
    glycan_pie_plot = glycan_pie_plot,
    glycan_bar_plot = glycan_bar_plot,
    enrichment_results = enrichment_results,
    glycan_data = glycan_data,
    significant_genes = significant_genes,
    unique_proteins = unique_proteins
  ))
}

# KEGG分析函数
kegg_analysis <- function(expr_df, group_df, glycan_type_df) {
  library(dplyr)
  library(tidyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(patchwork)
  library(gridExtra)
  # --- 数据预处理，与LDA/GO分析一致 ---
  colnames(expr_df) <- as.character(expr_df[1,])
  expr_df <- expr_df[-1,]
  gene_expression_long <- tidyr::gather(expr_df, key = "sample", value = "expression", -Glycopeptides)
  merged_data <- merge(gene_expression_long, group_df, by.x = "sample", by.y = "sample")
  wide_data <- tidyr::spread(merged_data, key = Glycopeptides, value = expression)
  wide_data[ , 3:ncol(wide_data)] <- lapply(wide_data[ , 3:ncol(wide_data)], as.numeric)
  # --- 计算均值比，筛选显著性基因 ---
  data_long <- wide_data %>%
    gather(key = "gene", value = "expression", -sample, -group)
  data_long$expression <- as.numeric(data_long$expression)
  group_means <- data_long %>%
    group_by(group, gene) %>%
    summarize(mean_expression = mean(expression, na.rm = TRUE), .groups = 'drop')
  group_means_wide <- group_means %>%
    spread(key = group, value = mean_expression)
  differences <- group_means_wide %>%
    mutate(CRLM_minus_CRC = CRLM / CRC,
           CRLM_minus_LIHC = CRLM / LIHC) %>%
    dplyr::select(gene, CRLM_minus_CRC, CRLM_minus_LIHC)
  regions_df <- differences %>%
    mutate(region = case_when(
      CRLM_minus_LIHC < 0.83 & CRLM_minus_CRC < 0.83 ~ "Bottom-Left",
      CRLM_minus_LIHC < 0.83 & CRLM_minus_CRC >= 0.83 & CRLM_minus_CRC <= 1.2 ~ "Left-Center",
      CRLM_minus_LIHC < 0.83 & CRLM_minus_CRC > 1.2 ~ "Top-Left",
      CRLM_minus_LIHC >= 0.83 & CRLM_minus_LIHC <= 1.2 & CRLM_minus_CRC < 0.83 ~ "Bottom-Center",
      CRLM_minus_LIHC >= 0.83 & CRLM_minus_LIHC <= 1.2 & CRLM_minus_CRC >= 0.83 & CRLM_minus_CRC <= 1.2 ~ "Center",
      CRLM_minus_LIHC >= 0.83 & CRLM_minus_LIHC <= 1.2 & CRLM_minus_CRC > 1.2 ~ "Top-Center",
      CRLM_minus_LIHC > 1.2 & CRLM_minus_CRC < 0.83 ~ "Bottom-Right",
      CRLM_minus_LIHC > 1.2 & CRLM_minus_CRC >= 0.83 & CRLM_minus_CRC <= 1.2 ~ "Right-Center",
      CRLM_minus_LIHC > 1.2 & CRLM_minus_CRC > 1.2 ~ "Top-Right"
    ))
  # 获取非中心区域的基因
  different_gene_mean_calculate_regions <- regions_df %>% 
    filter(region != "Center") %>% 
    pull(gene)
  # --- 热图分析获取combined_genes ---
  crlm_means <- data_long %>%
    filter(group == "CRLM") %>%
    group_by(gene) %>%
    summarize(mean_expression = mean(expression, na.rm = TRUE))
  crc_data <- wide_data %>% filter(group == "CRC")
  crc_data_long <- crc_data %>% gather(key = "gene", value = "expression", -sample, -group)
  crc_data_long$expression <- as.numeric(crc_data_long$expression)
  crc_vs_crlm <- crc_data_long %>%
    left_join(crlm_means, by = "gene") %>%
    mutate(expression_diff = mean_expression / expression) %>%
    dplyr::select(sample, gene, expression_diff)
  crc_vs_crlm_wide <- spread(crc_vs_crlm, key = gene, value = expression_diff)
  crc_expression_diff_matrix <- as.matrix(crc_vs_crlm_wide[,-1])
  rownames(crc_expression_diff_matrix) <- crc_vs_crlm_wide$sample
  crc_gene_filter <- apply(crc_expression_diff_matrix, 2, function(x) {
    sum(x > 1.2) >= length(x)/2 || sum(x < 0.83) >= length(x)/2
  })
  crc_filtered_genes <- colnames(crc_expression_diff_matrix)[crc_gene_filter]
  lihc_data <- wide_data %>% filter(group == "LIHC")
  lihc_data_long <- lihc_data %>% gather(key = "gene", value = "expression", -sample, -group)
  lihc_data_long$expression <- as.numeric(lihc_data_long$expression)
  lihc_vs_crlm <- lihc_data_long %>%
    left_join(crlm_means, by = "gene") %>%
    mutate(expression_diff = mean_expression / expression) %>%
    dplyr::select(sample, gene, expression_diff)
  lihc_vs_crlm_wide <- spread(lihc_vs_crlm, key = gene, value = expression_diff)
  lihc_expression_diff_matrix <- as.matrix(lihc_vs_crlm_wide[,-1])
  rownames(lihc_expression_diff_matrix) <- lihc_vs_crlm_wide$sample
  lihc_gene_filter <- apply(lihc_expression_diff_matrix, 2, function(x) {
    sum(x > 1.2) >= length(x)/2 || sum(x < 0.83) >= length(x)/2
  })
  lihc_filtered_genes <- colnames(lihc_expression_diff_matrix)[lihc_gene_filter]
  # 合并筛选基因
  combined_genes <- union(crc_filtered_genes, lihc_filtered_genes)
  # 统计检验获取显著性差异基因
  different_genes_from_both <- union(different_gene_mean_calculate_regions, combined_genes)
  significant_genes <- c()
  for (gene in different_genes_from_both) {
    gene_expression <- wide_data[, c("sample", "group", gene)]
    colnames(gene_expression) <- c("sample", "group", "expression")
    t_test_crc_crlm <- t.test(expression ~ group, data = subset(gene_expression, group %in% c("CRC", "CRLM")))
    t_test_crlm_lihc <- t.test(expression ~ group, data = subset(gene_expression, group %in% c("CRLM", "LIHC")))
    if (t_test_crc_crlm$p.value < 0.05 || t_test_crlm_lihc$p.value < 0.05) {
      significant_genes <- c(significant_genes, gene)
    }
  }
  significant_genes <- unique(significant_genes)
  # --- 提取蛋白ID并转换为ENTREZ ID ---
  sig_df <- as.data.frame(significant_genes)
  sig_df <- tidyr::separate(sig_df, significant_genes, into = c("peptides", "glycan", "protein"), sep = "-")
  unique_proteins <- unique(sig_df$protein)
  entrez_ids <- bitr(unique_proteins, fromType = "UNIPROT", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
  if (nrow(entrez_ids) == 0) {
    return(list(barplot = NULL, emapplot = NULL, cnetplot = NULL, kegg_table = NULL, kegg_obj = NULL, combined_plot = NULL))
  }
  # --- KEGG富集分析 ---
  kegg_enrichment <- tryCatch({
    enrichKEGG(gene = entrez_ids$ENTREZID,
               organism = 'hsa',
               pvalueCutoff = 0.05)
  }, error = function(e) NULL)
  if (is.null(kegg_enrichment) || nrow(kegg_enrichment@result) == 0) {
    return(list(barplot = NULL, emapplot = NULL, cnetplot = NULL, kegg_table = NULL, kegg_obj = NULL, combined_plot = NULL))
  }
  # 转换ENTREZ ID为基因符号
  if (nrow(kegg_enrichment@result) > 0) {
    kegg_enrichment@result$geneID <- sapply(kegg_enrichment@result$geneID, function(x) {
      gene_ids <- unlist(strsplit(x, "/"))
      gene_symbols <- entrez_ids$SYMBOL[match(gene_ids, entrez_ids$ENTREZID)]
      paste(gene_symbols, collapse = "/")
    })
    # 计算KEGG术语相似性矩阵
    kegg_enrichment <- pairwise_termsim(kegg_enrichment)
    # barplot
    barplot_obj <- barplot(kegg_enrichment, showCategory = 20, title = "KEGG Pathway Enrichment")
    # emapplot
    emapplot_obj <- tryCatch({ emapplot(kegg_enrichment, showCategory = 20) }, error = function(e) NULL)
    # cnetplot
    cnetplot_obj <- tryCatch({ cnetplot(kegg_enrichment, showCategory = 10) }, error = function(e) NULL)
    kegg_table <- as.data.frame(kegg_enrichment@result)
    
    # 生成组合图表
    if (!is.null(barplot_obj) && !is.null(emapplot_obj) && !is.null(cnetplot_obj)) {
      # 如果有三个图，创建2x2布局（第四个位置留空）
      combined_plot <- (barplot_obj + emapplot_obj) / (cnetplot_obj + plot_spacer()) +
        plot_layout(guides = "collect") +
        plot_annotation(
          title = "KEGG Pathway Analysis Results",
          theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
        )
    } else if (!is.null(barplot_obj) && !is.null(emapplot_obj)) {
      # 如果只有两个图，创建1x2布局
      combined_plot <- barplot_obj + emapplot_obj +
        plot_layout(guides = "collect") +
        plot_annotation(
          title = "KEGG Pathway Analysis Results",
          theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
        )
    } else if (!is.null(barplot_obj)) {
      # 如果只有一个图
      combined_plot <- barplot_obj +
        plot_annotation(
          title = "KEGG Pathway Analysis Results",
          theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
        )
    } else {
      combined_plot <- NULL
    }
  } else {
    barplot_obj <- NULL
    emapplot_obj <- NULL
    cnetplot_obj <- NULL
    kegg_table <- NULL
    combined_plot <- NULL
  }
  return(list(
    barplot = barplot_obj,
    emapplot = emapplot_obj,
    cnetplot = cnetplot_obj,
    kegg_table = kegg_table,
    kegg_obj = kegg_enrichment,
    combined_plot = combined_plot
  ))
}

# ROC曲线分析函数
roc_analysis <- function(expr_df, group_df, glycan_type_df) {
  library(dplyr)
  library(tidyr)
  library(pROC)
  library(caret)
  library(ggplot2)
  library(patchwork)
  # --- 数据预处理，与KEGG分析一致 ---
  colnames(expr_df) <- as.character(expr_df[1,])
  expr_df <- expr_df[-1,]
  gene_expression_long <- tidyr::gather(expr_df, key = "sample", value = "expression", -Glycopeptides)
  merged_data <- merge(gene_expression_long, group_df, by.x = "sample", by.y = "sample")
  wide_data <- tidyr::spread(merged_data, key = Glycopeptides, value = expression)
  wide_data[ , 3:ncol(wide_data)] <- lapply(wide_data[ , 3:ncol(wide_data)], as.numeric)
  # --- 计算均值比，筛选显著性基因 ---
  data_long <- wide_data %>%
    gather(key = "gene", value = "expression", -sample, -group)
  data_long$expression <- as.numeric(data_long$expression)
  group_means <- data_long %>%
    group_by(group, gene) %>%
    summarize(mean_expression = mean(expression, na.rm = TRUE), .groups = 'drop')
  group_means_wide <- group_means %>%
    spread(key = group, value = mean_expression)
  differences <- group_means_wide %>%
    mutate(CRLM_minus_CRC = CRLM / CRC,
           CRLM_minus_LIHC = CRLM / LIHC) %>%
    dplyr::select(gene, CRLM_minus_CRC, CRLM_minus_LIHC)
  regions_df <- differences %>%
    mutate(region = case_when(
      CRLM_minus_LIHC < 0.83 & CRLM_minus_CRC < 0.83 ~ "Bottom-Left",
      CRLM_minus_LIHC < 0.83 & CRLM_minus_CRC >= 0.83 & CRLM_minus_CRC <= 1.2 ~ "Left-Center",
      CRLM_minus_LIHC < 0.83 & CRLM_minus_CRC > 1.2 ~ "Top-Left",
      CRLM_minus_LIHC >= 0.83 & CRLM_minus_LIHC <= 1.2 & CRLM_minus_CRC < 0.83 ~ "Bottom-Center",
      CRLM_minus_LIHC >= 0.83 & CRLM_minus_LIHC <= 1.2 & CRLM_minus_CRC >= 0.83 & CRLM_minus_CRC <= 1.2 ~ "Center",
      CRLM_minus_LIHC >= 0.83 & CRLM_minus_LIHC <= 1.2 & CRLM_minus_CRC > 1.2 ~ "Top-Center",
      CRLM_minus_LIHC > 1.2 & CRLM_minus_CRC < 0.83 ~ "Bottom-Right",
      CRLM_minus_LIHC > 1.2 & CRLM_minus_CRC >= 0.83 & CRLM_minus_CRC <= 1.2 ~ "Right-Center",
      CRLM_minus_LIHC > 1.2 & CRLM_minus_CRC > 1.2 ~ "Top-Right"
    ))
  # 获取LIHC vs CRLM显著性基因
  crlm_means <- data_long %>%
    filter(group == "CRLM") %>%
    group_by(gene) %>%
    summarize(mean_expression = mean(expression, na.rm = TRUE))
  lihc_data <- wide_data %>% filter(group == "LIHC")
  lihc_data_long <- lihc_data %>% gather(key = "gene", value = "expression", -sample, -group)
  lihc_data_long$expression <- as.numeric(lihc_data_long$expression)
  lihc_vs_crlm <- lihc_data_long %>%
    left_join(crlm_means, by = "gene") %>%
    mutate(expression_diff = mean_expression / expression) %>%
    dplyr::select(sample, gene, expression_diff)
  lihc_vs_crlm_wide <- spread(lihc_vs_crlm, key = gene, value = expression_diff)
  lihc_expression_diff_matrix <- as.matrix(lihc_vs_crlm_wide[,-1])
  rownames(lihc_expression_diff_matrix) <- lihc_vs_crlm_wide$sample
  lihc_gene_filter <- apply(lihc_expression_diff_matrix, 2, function(x) {
    sum(x > 1.2) >= length(x)/2 || sum(x < 0.83) >= length(x)/2
  })
  final_significant_genes_lihc <- colnames(lihc_expression_diff_matrix)[lihc_gene_filter]
  # --- 构建ROC分析数据 ---
  data <- wide_data[, final_significant_genes_lihc, drop = FALSE]
  data$group <- wide_data$group
  data <- subset(data, data$group != "CRC")
  data <- data %>% mutate(group = ifelse(group == "LIHC", 0, ifelse(group == "CRLM", 1, group)))
  data$group <- as.numeric(data$group)
  group_vec <- data$group
  # 计算AUC并排序
  auc_results <- data.frame(Gene = character(), AUC = numeric(), stringsAsFactors = FALSE)
  for (gene in setdiff(names(data), "group")) {
    roc_curve <- roc(group_vec, data[[gene]])
    auc_value <- auc(roc_curve)
    auc_results <- rbind(auc_results, data.frame(Gene = gene, AUC = auc_value))
  }
  auc_results <- auc_results[order(auc_results$AUC, decreasing = TRUE), ]
  top_genes <- auc_results[1:min(5, nrow(auc_results)), "Gene"]
  # 联合模型
  selected_data <- data[, c(top_genes, "group")]
  model <- glm(group ~ ., data = selected_data, family = binomial)
  predicted_probabilities <- predict(model, type = "response")
  combined_roc_curve <- roc(selected_data$group, predicted_probabilities)
  # 联合ROC ggplot
  combined_roc_df <- data.frame(
    tpr = rev(combined_roc_curve$sensitivities),
    fpr = rev(1 - combined_roc_curve$specificities)
  )
  combined_roc_plot <- ggplot(combined_roc_df, aes(x = fpr, y = tpr)) +
    geom_line(color = "blue", size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray") +
    labs(title = "联合基因ROC曲线", x = "False positive rate", y = "True positive rate") +
    annotate("text", x = 0.8, y = 0.2, label = paste("AUC =", round(auc(combined_roc_curve), 2)), color = "blue") +
    theme_classic()
  # 单基因ROC ggplot
  single_roc_plots <- list()
  for (gene in top_genes) {
    single_roc_curve <- roc(group_vec, data[[gene]])
    single_roc_df <- data.frame(
      tpr = rev(single_roc_curve$sensitivities),
      fpr = rev(1 - single_roc_curve$specificities)
    )
    auc_val <- auc(single_roc_curve)
    p <- ggplot(single_roc_df, aes(x = fpr, y = tpr)) +
      geom_line(color = "blue", size = 1.2) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray") +
      labs(title = paste("ROC曲线", gene), x = "False positive rate", y = "True positive rate") +
      annotate("text", x = 0.8, y = 0.2, label = paste("AUC =", round(auc_val, 2)), color = "blue") +
      theme_classic()
    single_roc_plots[[gene]] <- p
  }
  # 合并单基因ROC为一张多图
  single_roc_patch <- wrap_plots(single_roc_plots, ncol = 2)
  return(list(
    combined_roc_plot = combined_roc_plot,
    single_roc_patch = single_roc_patch,
    auc_table = auc_results,
    top_genes = top_genes
  ))
}

# 火山图分析函数
volcano_analysis <- function(expr_df, group_df, glycan_type_df) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(patchwork)
  
  # --- 数据预处理，与其他分析一致 ---
  colnames(expr_df) <- as.character(expr_df[1,])
  expr_df <- expr_df[-1,]
  gene_expression_long <- tidyr::gather(expr_df, key = "sample", value = "expression", -Glycopeptides)
  merged_data <- merge(gene_expression_long, group_df, by.x = "sample", by.y = "sample")
  wide_data <- tidyr::spread(merged_data, key = Glycopeptides, value = expression)
  wide_data[, 3:ncol(wide_data)] <- lapply(wide_data[, 3:ncol(wide_data)], as.numeric)
  
  # 将 group 转换为因子类型
  wide_data$group <- factor(wide_data$group, levels = c("CRLM", "CRC", "LIHC"))
  
  # --- CRLM vs CRC 火山图 ---
  crlm_vs_crc <- wide_data %>%
    dplyr::select(-sample, -group) %>%
    bind_cols(group = wide_data$group) %>%
    gather(key = "gene", value = "expression", -group) %>%
    group_by(gene) %>%
    summarise(
      logFC = log2(mean(expression[group == "CRLM"]) / mean(expression[group == "CRC"])),
      p_value = t.test(expression[group == "CRLM"], expression[group == "CRC"])$p.value
    )
  
  crlm_vs_crc <- crlm_vs_crc %>%
    mutate(
      significance = case_when(
        p_value < 0.05 & abs(logFC) > 0.263034 ~ "Significant",
        TRUE ~ "Not Significant"
      )
    )
  
  # CRLM vs CRC 火山图1（标注所有显著基因）
  volcano_crc_1 <- ggplot(crlm_vs_crc, aes(x = logFC, y = -log10(p_value), color = significance)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
    theme_minimal() +
    ggtitle("Volcano Plot: CRLM vs CRC") +
    xlab("log2 Fold Change") +
    ylab("-log10(p-value)") +
    theme(legend.position = "top") +
    geom_text_repel(aes(label = ifelse(significance == "Significant", gene, "")), 
                    box.padding = 0.5, point.padding = 0.5, segment.size = 0.5, 
                    segment.color = "black", show.legend = FALSE)
  
  # CRLM vs CRC 火山图2（标注高表达显著基因）
  volcano_crc_2 <- ggplot(crlm_vs_crc, aes(x = logFC, y = -log10(p_value), color = significance)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
    theme_minimal() +
    ggtitle("Volcano Plot: CRLM vs CRC") +
    xlab("log2 Fold Change") +
    ylab("-log10(p-value)") +
    theme(legend.position = "top") +
    geom_text_repel(aes(label = ifelse(significance == "Significant" & abs(logFC) > 0.585, gene, "")), 
                    box.padding = 0.5, point.padding = 0.5, segment.size = 0.5, 
                    segment.color = "black", show.legend = FALSE)
  
  # --- CRLM vs LIHC 火山图 ---
  crlm_vs_lihc <- wide_data %>%
    dplyr::select(-sample, -group) %>%
    bind_cols(group = wide_data$group) %>%
    gather(key = "gene", value = "expression", -group) %>%
    group_by(gene) %>%
    summarise(
      logFC = log2(mean(expression[group == "CRLM"]) / mean(expression[group == "LIHC"])),
      p_value = t.test(expression[group == "CRLM"], expression[group == "LIHC"])$p.value
    )
  
  crlm_vs_lihc <- crlm_vs_lihc %>%
    mutate(
      significance = case_when(
        p_value < 0.05 & abs(logFC) > 0.263034 ~ "Significant",
        TRUE ~ "Not Significant"
      )
    )
  
  # CRLM vs LIHC 火山图1（标注所有显著基因）
  volcano_lihc_1 <- ggplot(crlm_vs_lihc, aes(x = logFC, y = -log10(p_value), color = significance)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
    theme_minimal() +
    ggtitle("Volcano Plot: CRLM vs LIHC") +
    xlab("log2 Fold Change") +
    ylab("-log10(p-value)") +
    theme(legend.position = "top") +
    geom_text_repel(aes(label = ifelse(significance == "Significant", gene, "")), 
                    box.padding = 0.5, point.padding = 0.5, segment.size = 0.5, 
                    segment.color = "black", show.legend = FALSE)
  
  # CRLM vs LIHC 火山图2（标注高表达显著基因）
  volcano_lihc_2 <- ggplot(crlm_vs_lihc, aes(x = logFC, y = -log10(p_value), color = significance)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
    theme_minimal() +
    ggtitle("Volcano Plot: CRLM vs LIHC") +
    xlab("log2 Fold Change") +
    ylab("-log10(p-value)") +
    theme(legend.position = "top") +
    geom_text_repel(aes(label = ifelse(significance == "Significant" & abs(logFC) > 0.585, gene, "")), 
                    box.padding = 0.5, point.padding = 0.5, segment.size = 0.5, 
                    segment.color = "black", show.legend = FALSE)
  
  # --- 蛋白质表达分析 ---
  # 定义要分析的蛋白质ID列表
  protein_ids <- c("P02763", "P01009", "P01008", "P05155", "P12109", "P02461", 
                   "Q99715", "P05997", "Q05707", "P12110", "P12111", "P06731", 
                   "P13688", "P07339")
  
  protein_plots <- list()
  
  for (protein_id in protein_ids) {
    # 提取包含该蛋白质的数据
    protein_data <- wide_data %>%
      dplyr::select(sample, group, contains(paste0("-", protein_id)))
    
    if (ncol(protein_data) > 2) {  # 确保有该蛋白质的数据
      protein_data <- protein_data %>%
        mutate(Sum = rowSums(dplyr::select(., -sample, -group)))
      
      comparisons <- list(c("CRLM", "LIHC"), c("CRLM", "CRC"))
      
      # 小提琴图
      violin_plot <- ggplot(protein_data, aes(x = group, y = Sum, fill = group)) +
        geom_violin(trim = FALSE, color = "black") +
        geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
        theme_classic() +
        labs(title = paste0("原验证数据-", protein_id, "_IGP"),
             x = "Group",
             y = "Expression Level",
             fill = "Group") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif")
      
      # 箱线图
      box_plot <- ggplot(protein_data, aes(x = group, y = Sum, fill = group)) +
        geom_boxplot() +
        geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
        theme_classic() +
        labs(title = paste0("原验证数据-", protein_id, "_IGP"),
             x = "Group",
             y = "Expression Level",
             fill = "Group") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif")
      
      protein_plots[[paste0(protein_id, "_violin")]] <- violin_plot
      protein_plots[[paste0(protein_id, "_box")]] <- box_plot
    }
  }
  
  # --- 聚类数据准备 ---
  df_long_for_cluster <- wide_data %>%
    pivot_longer(cols = -c(sample, group),
                 names_to = "variable",
                 values_to = "value")
  
  df_long_for_cluster <- df_long_for_cluster %>%
    group_by(group, variable) %>%
    summarize(mean_value = mean(value, na.rm = TRUE), .groups = 'drop')
  
  # 生成组合火山图
  volcano_combined_plot <- (volcano_crc_1 + volcano_crc_2) / (volcano_lihc_1 + volcano_lihc_2) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Volcano Plot Analysis Results",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
  
  return(list(
    # 火山图
    volcano_crc_1 = volcano_crc_1,
    volcano_crc_2 = volcano_crc_2,
    volcano_lihc_1 = volcano_lihc_1,
    volcano_lihc_2 = volcano_lihc_2,
    # 蛋白质表达图
    protein_plots = protein_plots,
    # 数据
    crlm_vs_crc_data = crlm_vs_crc,
    crlm_vs_lihc_data = crlm_vs_lihc,
    cluster_data = df_long_for_cluster,
    volcano_combined_plot = volcano_combined_plot
  ))
}

# 数据处理流程：批量加载、组内归一化、组间归一化
process_glycopeptide_files <- function(file_paths) {
  library(dplyr)
  library(tidyr)
  all_data_list <- list()
  # 1. 批量读取，若无raw.file则自动补充
  for (i in seq_along(file_paths)) {
    dat <- read.csv(file_paths[i], header = TRUE, sep = ",")
    # 检查raw.file列
    if (!("raw.file" %in% colnames(dat)) || all(is.na(dat$raw.file)) || all(dat$raw.file == "")) {
      dat$raw.file <- paste0("raw", i)
    }
    all_data_list[[i]] <- dat
  }
  # 2. 合并
  data_all <- do.call(rbind, all_data_list)
  # 3. 按raw.file分组做组内归一化
  list_all <- split(data_all, data_all$raw.file)
  zunei_normalization <- function(data) {
    data <- dplyr::filter(data, FDR <= 0.01) %>%
      dplyr::filter(NGLYCAN != "NA") %>%
      dplyr::filter(Protein.Name != "DECOY")
    data <- separate(data, Reporters, sep = ";", into = c("1","2","3","4","5","6","7","8","9","10"))
    data <- separate(data, "1", sep = "/", into = c("R_126","intensity_126_intensity"))
    data <- separate(data, "2", sep = "/", into = c("R_127C","intensity_127C_intensity"))
    data <- separate(data, "3", sep = "/", into = c("R_127N","intensity_127N_intensity"))
    data <- separate(data, "4", sep = "/", into = c("R_128C","intensity_128C_intensity"))
    data <- separate(data, "5", sep = "/", into = c("R_128N","intensity_128N_intensity"))
    data <- separate(data, "6", sep = "/", into = c("R_129C","intensity_129C_intensity"))
    data <- separate(data, "7", sep = "/", into = c("R_129N","intensity_129N_intensity"))
    data <- separate(data, "8", sep = "/", into = c("R_130C","intensity_130C_intensity"))
    data <- separate(data, "9", sep = "/", into = c("R_130N","intensity_130N_intensity"))
    data <- separate(data, "10", sep = "/", into = c("R_131","intensity_131_intensity"))
    data <- dplyr::select(data, -starts_with("R_"))
    data <- dplyr::select(data, Protein.Accession, Sequence, NGLYCAN, raw.file, starts_with("intensity_"))
    data <- separate(data, NGLYCAN, sep = ":", into = c("NGLYCAN", "NO"))
    data <- dplyr::select(data, -NO)
    data$Glycopeptides <- paste(data$Sequence, data$NGLYCAN, sep = "-")
    data$Protein.Accession <- sub(";.*", "", data$Protein.Accession)
    exp.names <- grep("^intensity", names(data), value = TRUE)
    data[exp.names] <- sapply(data[exp.names], as.numeric)
    data <- subset(data, data$intensity_131_intensity > 0)
    data[exp.names] <- sapply(data[exp.names], function(x) log2(x + 1))
    # 计算各组实验的相同糖肽的中值
    data_long <- pivot_longer(data, cols = starts_with("intensity_"), names_to = "sample", values_to = "before_scale_intensity")
    df_all <- split(data_long, data_long$sample)
    new_data <- data.frame()
    for (i in seq_along(df_all)) {
      glycopeptides_median <- data.frame()
      glycopeptides <- unique(df_all[[i]]$Glycopeptides)
      for (glyco in glycopeptides) {
        glyco_sub <- subset(df_all[[i]], df_all[[i]]$Glycopeptides == glyco)
        median_glyco_sub <- median(glyco_sub$before_scale_intensity, na.rm = TRUE)
        glyco_sub$before_scale_intensity <- rep(median_glyco_sub, length(glyco_sub$Glycopeptides))
        glycopeptides_median <- rbind(glycopeptides_median, glyco_sub)
      }
      new_data <- rbind(new_data, glycopeptides_median)
    }
    new_data <- unique(new_data)
    data_median_131_sub <- subset(new_data, new_data$sample == "intensity_131_intensity")
    data <- merge(new_data, data_median_131_sub, by = "Glycopeptides", all.x = TRUE)
    data$before_scale_intensity <- data$before_scale_intensity.x - data$before_scale_intensity.y
    # 恢复原始raw.file和sample信息
    data$raw.file <- data$raw.file.x
    data$sample <- data$sample.x
    data$Protein.Accession <- data$Protein.Accession.x
    data <- dplyr::select(data, Glycopeptides, Protein.Accession, raw.file, sample, before_scale_intensity)
    return(data)
  }
  for (z in seq_along(list_all)) {
    list_all[[z]] <- zunei_normalization(list_all[[z]])
  }
  all_data_after_zunei_normalization <- do.call(rbind, list_all)
  # 组间归一化
  channel_131_zujian <- subset(all_data_after_zunei_normalization, all_data_after_zunei_normalization$sample == "intensity_131_intensity")
  channel_131_zujian <- dplyr::select(channel_131_zujian, Glycopeptides, before_scale_intensity)
  colnames(channel_131_zujian) <- c("Glycopeptides", "intensity_median")
  glycopeptides <- unique(channel_131_zujian$Glycopeptides)
  new_131_glycopeptides <- data.frame()
  for (glyco in glycopeptides) {
    glyco_sub <- subset(channel_131_zujian, channel_131_zujian$Glycopeptides == glyco)
    median_glyco_sub <- median(glyco_sub$intensity_median)
    glyco_sub$intensity_median <- rep(median_glyco_sub, length(glyco_sub$Glycopeptides))
    new_131_glycopeptides <- rbind(new_131_glycopeptides, glyco_sub)
  }
  new_131_glycopeptides_zujian <- unique(new_131_glycopeptides)
  all_data_zujian_normalization <- merge(all_data_after_zunei_normalization, new_131_glycopeptides_zujian, by = "Glycopeptides", all.x = TRUE)
  all_data_zujian_normalization$final_scale_intensity <- all_data_zujian_normalization$before_scale_intensity + all_data_zujian_normalization$intensity_median
  # 只保留131通道以外的样本
  all_data_zujian_normalization <- subset(all_data_zujian_normalization, sample != "intensity_131_intensity")
  all_data_zujian_normalization$sample_id <- paste(all_data_zujian_normalization$raw.file, all_data_zujian_normalization$sample, sep = "_")
  all_data_zujian_normalization <- dplyr::select(all_data_zujian_normalization, Glycopeptides, Protein.Accession, raw.file, sample, final_scale_intensity, sample_id)
  return(list(
    normalized = all_data_zujian_normalization
  ))
}

# 严格复刻 FFPE_数据标准化_糖组分析_糖肽定性分析_250215.R 的数据处理流程
process_glycopeptide_files_strict <- function(file_paths, raw_file_names = NULL) {
  library(dplyr)
  library(tidyr)
  all_data_list <- list()
  n_files <- length(file_paths)
  if (is.null(raw_file_names)) {
    raw_file_names <- paste0("raw", seq_len(n_files))
  }
  for (i in seq_along(file_paths)) {
    dat <- read.csv(file_paths[i], header = TRUE, sep = ",")
    dat$raw.file <- rep(raw_file_names[i], nrow(dat))
    all_data_list[[i]] <- dat
  }
  data_all <- do.call(rbind, all_data_list)
  list_all <- split(data_all, data_all$raw.file)
  zunei_normalization <- function(data) {
    data <- dplyr::filter(data, FDR <= 0.01) %>%
      dplyr::filter(NGLYCAN != "NA") %>%
      dplyr::filter(Protein.Name != "DECOY")
    data <- separate(data, Reporters, sep = ";", into = c("1","2","3","4","5","6","7","8","9","10"))
    data <- separate(data, "1", sep = "/", into = c("R_126","intensity_126_intensity"))
    data <- separate(data, "2", sep = "/", into = c("R_127C","intensity_127C_intensity"))
    data <- separate(data, "3", sep = "/", into = c("R_127N","intensity_127N_intensity"))
    data <- separate(data, "4", sep = "/", into = c("R_128C","intensity_128C_intensity"))
    data <- separate(data, "5", sep = "/", into = c("R_128N","intensity_128N_intensity"))
    data <- separate(data, "6", sep = "/", into = c("R_129C","intensity_129C_intensity"))
    data <- separate(data, "7", sep = "/", into = c("R_129N","intensity_129N_intensity"))
    data <- separate(data, "8", sep = "/", into = c("R_130C","intensity_130C_intensity"))
    data <- separate(data, "9", sep = "/", into = c("R_130N","intensity_130N_intensity"))
    data <- separate(data, "10", sep = "/", into = c("R_131","intensity_131_intensity"))
    data <- dplyr::select(data, -starts_with("R_"))
    data <- dplyr::select(data, Protein.Accession, Sequence, NGLYCAN, raw.file, starts_with("intensity_"))
    data <- separate(data, NGLYCAN, sep = ":", into = c("NGLYCAN", "NO"))
    data <- dplyr::select(data, -NO)
    data$Glycopeptides <- paste(data$Sequence, data$NGLYCAN, sep = "-")
    data$Protein.Accession <- sub(";.*", "", data$Protein.Accession)
    exp.names <- grep("^intensity", names(data), value = TRUE)
    data[exp.names] <- sapply(data[exp.names], as.numeric)
    data <- subset(data, data$intensity_131_intensity > 0)
    data[exp.names] <- sapply(data[exp.names], function(x) log2(x + 1))
    data_long <- pivot_longer(data, cols = starts_with("intensity_"), names_to = "sample", values_to = "before_scale_intensity")
    df_all <- split(data_long, data_long$sample)
    new_data <- data.frame()
    for (i in seq_along(df_all)) {
      glycopeptides_median <- data.frame()
      glycopeptides <- unique(df_all[[i]]$Glycopeptides)
      for (glyco in glycopeptides) {
        glyco_sub <- subset(df_all[[i]], df_all[[i]]$Glycopeptides == glyco)
        median_glyco_sub <- median(glyco_sub$before_scale_intensity, na.rm = TRUE)
        glyco_sub$before_scale_intensity <- rep(median_glyco_sub, length(glyco_sub$Glycopeptides))
        glycopeptides_median <- rbind(glycopeptides_median, glyco_sub)
      }
      new_data <- rbind(new_data, glycopeptides_median)
    }
    new_data <- unique(new_data)
    data_median_131_sub <- subset(new_data, new_data$sample == "intensity_131_intensity")
    data <- merge(new_data, data_median_131_sub, by = "Glycopeptides", all.x = TRUE)
    data$before_scale_intensity <- data$before_scale_intensity.x - data$before_scale_intensity.y
    return(data)
  }
  for (z in seq_along(list_all)) {
    list_all[[z]] <- zunei_normalization(list_all[[z]])
  }
  all_data_after_zunei_normalization <- do.call(rbind, list_all)
  channel_131_zujian <- subset(all_data_after_zunei_normalization, all_data_after_zunei_normalization$sample.x == "intensity_131_intensity")
  channel_131_zujian <- dplyr::select(channel_131_zujian, Glycopeptides, before_scale_intensity.x)
  colnames(channel_131_zujian) <- c("Glycopeptides", "intensity_median")
  glycopeptides <- unique(channel_131_zujian$Glycopeptides)
  new_131_glycopeptides <- data.frame()
  for (glyco in glycopeptides) {
    glyco_sub <- subset(channel_131_zujian, channel_131_zujian$Glycopeptides == glyco)
    median_glyco_sub <- median(glyco_sub$intensity_median)
    glyco_sub$intensity_median <- rep(median_glyco_sub, length(glyco_sub$Glycopeptides))
    new_131_glycopeptides <- rbind(new_131_glycopeptides, glyco_sub)
  }
  new_131_glycopeptides_zujian <- unique(new_131_glycopeptides)
  all_data_zujian_normalization <- merge(all_data_after_zunei_normalization, new_131_glycopeptides_zujian, by = "Glycopeptides", all.x = TRUE)
  all_data_zujian_normalization$final_scale_intensity <- all_data_zujian_normalization$before_scale_intensity + all_data_zujian_normalization$intensity_median
  all_data_zujian_normalization <- dplyr::select(all_data_zujian_normalization, Glycopeptides, Protein.Accession.x, raw.file.x, sample.x, final_scale_intensity)
  colnames(all_data_zujian_normalization) <- c("Glycopeptides", "Protein.Accession", "raw.file", "sample", "final_scale_intensity")
  all_data_zujian_normalization$sample_id <- paste(all_data_zujian_normalization$raw.file, all_data_zujian_normalization$sample, sep = "_")
  return(list(
    normalized = all_data_zujian_normalization
  ))
}