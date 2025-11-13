gc()
rm(list = ls())
library(dplyr)
library(tidyr)
setwd("D:/BiologicalCodeinR/proteomeExpert-proteomeExpert/FFPE_glycopeptides/FFPE_glycopeptides")

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






library(tidyr)
data_all=rbind(data2_1,data2_2,data2_3,data3_1,data3_2,data3_3,
               data4_1,data4_2,data4_3,data5_1,data5_2,data5_3,
               data6_1,data6_2,data6_3,data7_1,data7_2,data7_3,
               data8_1,data8_2,data8_3,data9_1,data9_2,data9_3)
list_all=split(data_all,data_all$raw.file)
colnames(list_all[[1]])

zunei_normalization=function(data){
  library(tidyr)
  ####初步筛选
  data=dplyr::filter(data, FDR<= 0.01)%>%
    dplyr::filter(NGLYCAN != "NA")%>%
    dplyr::filter(Protein.Name!="DECOY")
  ###先将实验组分开
  data=separate(data,Reporters,sep=";",into = c("1","2","3","4","5","6","7","8","9","10"))
  data=separate(data,"1",sep="/",into = c("R_126","intensity_126_intensity"))
  data=separate(data,"2",sep="/",into = c("R_127C","intensity_127C_intensity"))
  data=separate(data,"3",sep="/",into = c("R_127N","intensity_127N_intensity"))
  data=separate(data,"4",sep="/",into = c("R_128C","intensity_128C_intensity"))
  data=separate(data,"5",sep="/",into = c("R_128N","intensity_128N_intensity"))
  data=separate(data,"6",sep="/",into = c("R_129C","intensity_129C_intensity"))
  data=separate(data,"7",sep="/",into = c("R_129N","intensity_129N_intensity"))
  data=separate(data,"8",sep="/",into = c("R_130C","intensity_130C_intensity"))
  data=separate(data,"9",sep="/",into = c("R_130N","intensity_130N_intensity"))
  data=separate(data,"10",sep="/",into = c("R_131","intensity_131_intensity"))
  data=dplyr::select(data,-starts_with("R_"))
  data=dplyr::select(data,Protein.Accession,Sequence,NGLYCAN,raw.file, starts_with("intensity_"))
  data=separate(data,NGLYCAN,sep=":",into = c("NGLYCAN","NO"))
  data=dplyr::select(data,-NO)
  # 使用paste替代unite函数，避免版本兼容性问题
  data$Glycopeptides <- paste(data$Sequence, data$NGLYCAN, sep = "-")
  data$Protein.Accession = sub(";.*", "", data$Protein.Accession)
  exp.names = grep("^intensity", names(data), value = TRUE)
  data[exp.names] = sapply(data[exp.names], as.numeric)
  data=subset(data, data$intensity_131_intensity>0)
  range(data$intensity_131_intensity)
  exp.names = grep("^intensity", names(data), value = TRUE)
  data[exp.names] = sapply(data[exp.names], as.numeric)
  data[exp.names] = sapply(data[exp.names],  function(x) log2(x + 1))
  
  #####上面做完log2转换之后，计算各组实验的相同糖肽的中值，
  data_long=pivot_longer(data, cols = starts_with("intensity_"), names_to = "sample", values_to = "before_scale_intensity")
  df_all=split(data_long,data_long$sample)
  new_data=data.frame()
  lengthi=as.numeric(length(names(df_all)))
  i=1
  for (i in 1:lengthi){
    glycopeptides_median=data.frame()
    glycopeptides=unique(df_all[[i]]$Glycopeptides)
    lengthj=as.numeric(length(glycopeptides))
    j=1
    for (j in 1:lengthj){
      glyco_sub=subset(df_all[[i]],df_all[[i]]$Glycopeptides==glycopeptides[j])
      median_glyco_sub=median(glyco_sub$before_scale_intensity,na.rm = T)
      glyco_sub$before_scale_intensity=rep(median_glyco_sub,as.numeric(length(glyco_sub$Glycopeptides)))
      glycopeptides_median=rbind(glycopeptides_median,glyco_sub)
      j=j+1
    }
    new_data=rbind(new_data,glycopeptides_median)
    i=i+1
  }
  new_data=unique(new_data)
  data_median_131_sub=subset(new_data,new_data$sample=="intensity_131_intensity")
  data=merge(new_data,data_median_131_sub,by="Glycopeptides",all.x = T)
  data$before_scale_intensity=data$before_scale_intensity.x - data$before_scale_intensity.y
  return(data)
}


lengthz=as.numeric(length(names(list_all)))
z=1 
for (z in 1:lengthz){
  list_all[[z]]=zunei_normalization(list_all[[z]])
  z=z+1
}


i=1

all_data_after_zunei_normalization=data.frame()
for (i in 1:lengthz){
  all_data_after_zunei_normalization=rbind(all_data_after_zunei_normalization,list_all[[i]])
  i=i+1
}
unique(all_data_after_zunei_normalization$sample.x)
unique(all_data_after_zunei_normalization$raw.file.x)

####组间归一化
#####思路为：上述拿到的131数据取出，每个糖肽对应的数值改为组间中值
#####然后merge样品数据和131通道数据做加法
channel_131_zujian=subset(all_data_after_zunei_normalization,all_data_after_zunei_normalization$sample.x=="intensity_131_intensity")
channel_131_zujian=dplyr::select(channel_131_zujian,Glycopeptides,before_scale_intensity.x)
colnames(channel_131_zujian)=c("Glycopeptides","intensity_median")
range(channel_131_zujian$intensity_median)

result <-channel_131_zujian%>%
  group_by(Glycopeptides) %>%
  summarise(median_intensity = median(intensity_median, na.rm = TRUE))

write.csv(result,file = "新数据131_median.csv")



glycopeptides=unique(channel_131_zujian$Glycopeptides)
lengthi=as.numeric(length(glycopeptides))
new_131_glycopeptides=data.frame()
i=1
for (i in 1:lengthi){
  glyco_sub=subset(channel_131_zujian,channel_131_zujian$Glycopeptides==glycopeptides[i])
  median_glyco_sub=median(glyco_sub$intensity_median)
  glyco_sub$intensity_median=rep(median_glyco_sub,as.numeric(length(glyco_sub$Glycopeptides)))
  new_131_glycopeptides=rbind(new_131_glycopeptides,glyco_sub)
  i=i+1
}
new_131_glycopeptides_zujian=unique(new_131_glycopeptides)



all_data_zujian_normalization=merge(all_data_after_zunei_normalization,new_131_glycopeptides_zujian,by="Glycopeptides",all.x = T)
all_data_zujian_normalization$final_scale_intensity=all_data_zujian_normalization$before_scale_intensity+ all_data_zujian_normalization$intensity_median  ############组间归一化
colnames(all_data_zujian_normalization)
all_data_zujian_normalization=dplyr::select(all_data_zujian_normalization,Glycopeptides,Protein.Accession.x,raw.file.x,sample.x,final_scale_intensity)
colnames(all_data_zujian_normalization)=c("Glycopeptides","Protein.Accession","raw.file","sample","final_scale_intensity")


all_data_zujian_normalization$sample_id=paste(all_data_zujian_normalization$raw.file,all_data_zujian_normalization$sample,sep = "_")
write.csv(all_data_zujian_normalization,"all_data_zujian_normal.csv", row.names = FALSE)  # 添加row.names = FALSE







###################################################################################################
gc()
rm(list = ls())
library(dplyr)
library(tidyr)
setwd("D:/BiologicalCodeinR/proteomeExpert-proteomeExpert/FFPE_glycopeptides/FFPE_glycopeptides")
data=read.csv("all_data_zujian_normal.csv",header = T)
data$final_scale_intensity=2^data$final_scale_intensity
library(ggplot2)
data = select(data, -raw.file, -sample)  
colnames(data)
# 加载必要的包
library(tidyr)
library(dplyr)
length(unique(data$Glycopeptides))

#####################################################################################################
##糖组数据分析
data_for_glycomics = select(data, -Protein.Accession)  
print(data_for_glycomics)
data_for_glycomics=separate(data_for_glycomics,Glycopeptides,sep = "-",into = c("Peptides","Glycan"))
data_for_glycomics = select(data_for_glycomics, -Peptides)  
print(data_for_glycomics)
data_for_glycomics <- data_for_glycomics %>%
  group_by(Glycan, sample_id) %>%
  summarise(total_intensity = sum(final_scale_intensity))
print(data_for_glycomics)
patients=read.csv("FFPE_patients_reflect.csv",header = T)
print(patients)
print(data_for_glycomics)
data_for_glycomics=merge(data_for_glycomics,patients,by.x="sample_id",by.y="Sample")
data_for_glycomics = select(data_for_glycomics, -sample_id) 
colnames(data_for_glycomics)
data_for_glycomics = select(data_for_glycomics, -patients, -repeat.) 
df_wide_for_glycomics <- data_for_glycomics %>%
  pivot_wider(names_from = patients_repeat, values_from = total_intensity)
###计算有效值
# 将NA替换为数值NA，以便后续计算
df_wide_for_glycomics[df_wide_for_glycomics == "NA"] <- NA
# 计算每个基因在所有样品中的有效值比例，并合并回原始数据框
df_wide_for_glycomics <- df_wide_for_glycomics %>%
  rowwise() %>%
  mutate(valid_ratio = mean(as.numeric(!is.na(c_across(-Glycan)))))
numeric_values <- df_wide_for_glycomics$valid_ratio[is.numeric(df_wide_for_glycomics$valid_ratio)]
summary(df_wide_for_glycomics$valid_ratio)
# 保留有效值大于   的基因
df_filtered <- df_wide_for_glycomics %>%
  filter(valid_ratio > 0.5) %>%
  dplyr::select(-valid_ratio)
library(VIM)  # 用于KNN填补缺失值
#install.packages("VIM")
# 使用KNN填补缺失值
df_filled <- kNN(df_filtered, 
                 variable = setdiff(colnames(df_filtered), "Glycan")) %>%
             dplyr::select(Glycan, everything())  
colnames(df_filled)
df_filled=df_filled %>%
  dplyr::select(-matches("_imp$"))
colnames(df_filled)

df_long_filled_glycomics <- df_filled %>%
  pivot_longer(cols = -Glycan,  # Keep the 'Glycan' column fixed
               names_to = "Sample",  # New column to hold the names of the other columns
               values_to = "Value")  # New column to hold the values of the other columns

colnames(df_long_filled_glycomics)
df_long_filled_glycomics=separate(df_long_filled_glycomics,Sample,sep = "\\.",into = c("patients","cancertype"))
df_long_filled_glycomics=separate(df_long_filled_glycomics,cancertype,sep = "_",into = c("cancertype","repeat"))
calculate_molar_mass <- function(glyco_str) {
  # 提取糖基信息
  matches <- regmatches(glyco_str, regexec("N(\\d+)H(\\d+)F(\\d+)S(\\d+)G(\\d+)", glyco_str))
  
  # 获取每个糖基的数量
  N <- as.integer(matches[[1]][2])
  H <- as.integer(matches[[1]][3])
  F <- as.integer(matches[[1]][4])
  S <- as.integer(matches[[1]][5])
  G <- as.integer(matches[[1]][6])
  
  # 各个糖的摩尔质量（单位：g/mol）
  molar_mass_GlcNAc <- 221.24   # GlcNAc
  molar_mass_Mannose <- 180.16  # Mannose
  molar_mass_Fucose <- 164.16   # Fucose
  molar_mass_Sialic_Acid <- 309.24  # Sialic acid
  molar_mass_Glucose <- 180.16  # Glucose
  
  # 计算总质量
  total_mass <- (N * molar_mass_GlcNAc) + (H * molar_mass_Mannose) + (F * molar_mass_Fucose) +
    (S * molar_mass_Sialic_Acid) + (G * molar_mass_Glucose)
  
  # 计算糖苷键数量
  # 假设每种糖基之间都形成糖苷键，且一个糖基之间的键数为糖基数量减1（例如2个糖基形成1个键）
  glycosidic_bonds <- (N + H  + F + S+G - 1)
  
  # 每个糖苷键会失去一个水分子（18.015 g/mol）
  water_loss <- glycosidic_bonds * 18.015
  
  # 计算最终的摩尔质量，减去水分子的质量
  final_mass <- total_mass - water_loss
  
  return(final_mass)
}

# 在df_long_summary中应用这个计算函数，生成新的列
df_long_filled_glycomics <- df_long_filled_glycomics %>%
  mutate(Glyco_Molar_Mass = sapply(Glycan, calculate_molar_mass))
colnames(df_long_filled_glycomics)

###############存在显著性差异的
# 根据癌症类型筛选数据
df_crc <- df_long_filled_glycomics %>% filter(cancertype == "CRC")
df_crlm <- df_long_filled_glycomics %>% filter(cancertype == "CRLM")
df_lihc <- df_long_filled_glycomics %>% filter(cancertype == "LIHC")

# 定义函数来执行 t 检验并提取结果
perform_t_test <- function(df1, df2, glycan, glyco_mass, comparison) {
  # 执行 t 检验
  t_test_result <- t.test(df1$Value, df2$Value)
  
  # 返回包含 Glycan、Glyco_Molar_Mass、p 值和比较信息的数据框
  return(data.frame(
    Glycan = glycan,
    Glyco_Molar_Mass = glyco_mass,
    p_value = t_test_result$p.value,
    comparison = comparison
  ))
}

# 针对每种 Glycan 和 Glyco_Molar_Mass 执行 t 检验
results_crc_crlm <- lapply(unique(df_long_filled_glycomics$Glycan), function(glycan) {
  glyco_mass <- unique(df_long_filled_glycomics$Glyco_Molar_Mass[df_long_filled_glycomics$Glycan == glycan])
  perform_t_test(df_crc[df_crc$Glycan == glycan, ], df_crlm[df_crlm$Glycan == glycan, ], glycan, glyco_mass, "CRC vs CRLM")
})

results_crlm_lihc <- lapply(unique(df_long_filled_glycomics$Glycan), function(glycan) {
  glyco_mass <- unique(df_long_filled_glycomics$Glyco_Molar_Mass[df_long_filled_glycomics$Glycan == glycan])
  perform_t_test(df_crlm[df_crlm$Glycan == glycan, ], df_lihc[df_lihc$Glycan == glycan, ], glycan, glyco_mass, "CRLM vs LIHC")
})

# 合并结果
results_df <- bind_rows(results_crc_crlm, results_crlm_lihc)

# 筛选显著性差异（p 值小于 0.05）
significant_results <- results_df %>% filter(p_value < 0.05)

# 查看显著性结果
head(significant_results)
df_long_filled_glycomics <-df_long_filled_glycomics%>%
  group_by(Glyco_Molar_Mass, cancertype) %>%
  summarise(mean_intensity = mean(Value, na.rm = TRUE))

data_for_normalized_glycomics <- df_long_filled_glycomics %>%
  group_by(cancertype) %>%
  mutate(
    Mean_Value_normalized = (mean_intensity / max(mean_intensity)) * 100  # 按最高比例归一化
  ) %>%
  ungroup()
colnames(data_for_normalized_glycomics)
# 获取所有 Glyco_Molar_Mass 的唯一值并排序
x_ticks <- sort(unique(data_for_normalized_glycomics$Glyco_Molar_Mass))
colnames(data_for_normalized_glycomics)
# 绘制图表
ggplot(data_for_normalized_glycomics, aes(x = Glyco_Molar_Mass, y = Mean_Value_normalized)) +
  geom_segment(aes(xend = Glyco_Molar_Mass, yend = 0), size = 1) +  # 绘制竖线
  facet_wrap(~ cancertype, scales = "fixed", ncol = 1) +  # 按 CancerType 分面竖着排列，坐标轴固定
  scale_x_continuous(
    limits = c(570, 6500),  # 设置 x 轴范围
    breaks = x_ticks  # 使用 Glyco_Molar_Mass 的所有唯一值作为 x 轴刻度
  ) +
  #scale_y_continuous(limits = c(0, 100)) +  # 设置 y 轴范围为 25 到 100
  labs(x = "Glyco_Molar_Mass", y = "Mean_Value") +  # 设置标签
  theme_minimal() +  # 使用简洁的主题
  theme(
    strip.text = element_text(size = 12),  # 调整分面标签大小
    axis.text = element_text(size = 10),   # 调整坐标轴标签大小
    axis.title = element_text(size = 12)   # 调整坐标轴标题大小
  )

significant_glycos <- significant_results$Glyco_Molar_Mass
highlighted_ticks <- significant_glycos

ggplot(data_for_normalized_glycomics, aes(x = Glyco_Molar_Mass, y = Mean_Value_normalized)) +
  geom_segment(aes(xend = Glyco_Molar_Mass, yend = 0), size = 1) +  # 绘制竖线
  facet_wrap(~ cancertype, scales = "fixed", ncol = 1) +  # 按 CancerType 分面竖着排列，坐标轴固定
  scale_x_continuous(
    limits = c(570, 5000),  # 设置 x 轴范围
    #breaks = x_ticks  # 使用 Glyco_Molar_Mass 的所有唯一值作为 x 轴刻度
  ) +
  # scale_y_continuous(limits = c(25, 100)) +  # 设置 y 轴范围为 25 到 100
  labs(x = "Glyco_Molar_Mass", y = "Mean_Value") +  # 设置标签
  theme_minimal() +  # 使用简洁的主题
  theme(
    strip.text = element_text(size = 12),  # 调整分面标签大小
    axis.text = element_text(size = 10),   # 调整坐标轴标签大小
    axis.title = element_text(size = 12)   # 调整坐标轴标题大小
  ) +
  # 标注红色的 Glyco_Molar_Mass 值
  geom_text(data = subset(data_for_normalized_glycomics, Glyco_Molar_Mass %in% highlighted_ticks), 
            aes(label = Glyco_Molar_Mass), color = "red", size = 3, vjust = -0.5)




####################################################################################################
##糖肽定性分析
# 将数据框转换为宽格式
df_wide <- data %>%
  pivot_wider(names_from = sample_id, values_from = final_scale_intensity)
df_wide <- df_wide[, !grepl("131", names(df_wide))]
# 查看结果
print(df_wide)
unique(data$sample_id)
colnames(df_wide)
df_wide$Glycopeptides=paste(df_wide$Glycopeptides,df_wide$Protein.Accession, sep = "-")
df_wide=df_wide %>% select(-Protein.Accession)
patients=read.csv("FFPE_patients_reflect.csv",header = T)
###列名替换成病人信息
# 提取旧的列名
old_colnames <- colnames(df_wide)
# 提取新的列名和样本名
samples <- patients$Sample
new_colnames <- patients$patients_repeat
# 创建一个命名向量用于替换列名，保持列名顺序
name_map <- setNames(new_colnames, samples)
# 替换 df_wide 中匹配到的列名
for (i in 1:length(old_colnames)) {
  if (old_colnames[i] %in% names(name_map)) {
    old_colnames[i] <- name_map[[old_colnames[i]]]
  }
}
# 应用新的列名到 df_wide
colnames(df_wide) <- old_colnames

# 将宽数据格式转换为长数据格式
df_long <- df_wide %>%
  pivot_longer(
    cols = -c(Glycopeptides),  # 指定要转换的列
    names_to = "Sample",      # 新列的名称，用于存储旧列的名称
    values_to = "Value"       # 新列的名称，用于存储旧列的值
  )
df_long=df_long%>%filter(!is.na(Value))

# 绘制抖动散点图加箱线图
p=ggplot(df_long, aes(x = Sample, y = Value)) +
  geom_boxplot() +  # 绘制箱线图，隐藏异常值
 # geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +  # 绘制抖动散点图
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # 调整x轴标签角度
  labs(title = "组间归一化", x = "Sample", y = "Value")+theme_classic()
#ggsave(p,file="组间归一化.pdf",width=25,height=8)
print(p)
########################################################################################################
df_for_dingxing_analysis=df_wide
glycan_type=read.csv("glycan-glycantype.csv",header = T)
df_for_dingxing_analysis$glycopeptides=df_for_dingxing_analysis$Glycopeptides
df_for_dingxing_analysis <- df_for_dingxing_analysis %>%
  separate(Glycopeptides, into = c("peptides", "glycan", "Protein"), sep = "-", extra = "merge")
df_for_dingxing_analysis=merge(df_for_dingxing_analysis,glycan_type)
last_col <- tail(names(df_for_dingxing_analysis), n = 1)

# 获取前三列的列名
first_three_cols <- c("glycan", "peptides", "Protein")

# 获取中间列的列名，如果存在中间列
if (ncol(df_for_dingxing_analysis) > 4) {
  middle_cols <- names(df_for_dingxing_analysis)[4:(ncol(df_for_dingxing_analysis)-1)]
} else {
  middle_cols <- character(0)
}
# 重新排列列顺序，将最后一列移动到第四列
new_col_order <- c(first_three_cols, last_col, middle_cols)

# 重排列顺序
df_for_dingxing_analysis <- df_for_dingxing_analysis[, new_col_order]


# 加载必要的包
library(dplyr)
library(tidyr)

df_for_dingxing_analysis_long <- pivot_longer(df_for_dingxing_analysis, 
                                             cols = -c(glycan, peptides, Protein, glycan.type, glycopeptides),
                                             names_to = "patient_info", 
                                             values_to = "value")
df_for_dingxing_analysis_long <- subset(df_for_dingxing_analysis_long, !is.na(df_for_dingxing_analysis_long$value))
unique(df_for_dingxing_analysis_long$patient_info)
df_for_dingxing_analysis_long <- df_for_dingxing_analysis_long %>%
  separate(patient_info, into = c("patients", "repeat"), sep = "_", extra = "merge")
df_for_dingxing_analysis_long <- df_for_dingxing_analysis_long %>%
  separate(patients, into = c("patients", "cancer_type"), sep = "\\.", extra = "merge")

# 加载必要的包
library(dplyr)
library(ggplot2)
unique(df_for_dingxing_analysis_long$glycan.type)
unique(df_for_dingxing_analysis_long$cancer_type)
# 假设数据框已经加载为 df_for_dingxing_analysis_long
table(df_for_dingxing_analysis_long$cancer_type)

# 按癌症类型分组并计算 unique 数目
unique_counts_by_cancer_type <- df_for_dingxing_analysis_long %>%
  filter(value > 0) %>%  # 过滤掉 value 小于等于 0 的数据
  group_by(cancer_type) %>%
  summarise(
    unique_proteins = n_distinct(Protein),
    unique_glycans = n_distinct(glycan),
    unique_peptides = n_distinct(peptides),
    unique_glycopeptides = n_distinct(glycopeptides)
  ) %>%
  pivot_longer(cols = c(unique_proteins, unique_glycans, unique_peptides, unique_glycopeptides), 
               names_to = "unique_type", values_to = "count")


# 查看统计结果
print(unique_counts_by_cancer_type)

# 使用 ggplot2 进行可视化
ggplot(unique_counts_by_cancer_type, aes(x = cancer_type, y = count, fill = unique_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Unique Counts by Cancer Type",
       x = "Cancer Type",
       y = "Count",
       fill = "Unique Type") +
  geom_text(aes(label = count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3) +
  theme_classic()

# 计算每个肽段上修饰的glycan类别数目
peptide_glycan=df_for_dingxing_analysis_long[,c("glycan", "peptides", "glycan.type")]
peptide_glycan=unique(peptide_glycan)
glycan_counts_per_peptide <- peptide_glycan %>%
  group_by(peptides) %>%
  summarise(unique_glycans = n_distinct(glycan))

# 查看统计结果
print(glycan_counts_per_peptide)

# 如果需要将结果保存为数据框
glycan_counts_per_peptide_df <- as.data.frame(glycan_counts_per_peptide)

# 使用 ggplot2 进行可视化 - 直方图
ggplot(glycan_counts_per_peptide_df, aes(x = unique_glycans)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Unique Glycan Types per Peptide",
       x = "Number of Unique Glycan Types",
       y = "Frequency") +
  theme_classic()

# 使用 ggplot2 进行可视化 - 密度图
ggplot(glycan_counts_per_peptide_df, aes(x = unique_glycans)) +
  geom_density(fill = "blue", alpha = 0.7) +
  labs(title = "Density of Unique Glycan Types per Peptide",
       x = "Number of Unique Glycan Types",
       y = "Density") +
  theme_classic()


# 计算糖肽类型
glycopeptide_glycan=df_for_dingxing_analysis_long[,c("glycan", "peptides", "glycan.type")]
glycopeptide_glycan$glycopep=paste(glycopeptide_glycan$peptides,glycopeptide_glycan$glycan,sep = "_")
glycopeptide_glycan=unique(glycopeptide_glycan)
# 计算每个 glycan.type 对应的 glycopeptides 数目
glycan_type_counts <- glycopeptide_glycan %>%
  group_by(glycan.type) %>%
  summarise(glycopeptides_count = n())
# 查看统计结果
print(glycan_type_counts)

# 如果需要将结果保存为数据框
glycan_type_counts_df <- as.data.frame(glycan_type_counts)

ggplot(glycan_type_counts_df, aes(x = glycan.type, y = glycopeptides_count)) +
  geom_bar(stat = "identity", fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Number of Glycopeptides per Glycan Type",
       x = "Glycan Type",
       y = "Number of Glycopeptides") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 旋转 x 轴标签使其可读


glyco_glycan_sites=merge(glycopeptide_glycan,glycan_counts_per_peptide_df)
colnames(glyco_glycan_sites)
colnames(glyco_glycan_sites)=c("glycan","peptides","glycan.type","glycopep","sites")
unique(glyco_glycan_sites$glycan.type)
glyco_glycan_sites=glyco_glycan_sites %>% select(-glycopep)
colnames(glyco_glycan_sites)
colnames(glyco_glycan_sites)=c("peptides","glycan","glycan.type","sites")
#####cancer glycantype sites桑葚图分析
df_for_dingxing_analysis_long_sangshen=merge(df_for_dingxing_analysis_long,glyco_glycan_sites)
table(df_for_dingxing_analysis_long_sangshen$sites)
# 按列名删除指定列
df_for_dingxing_analysis_long_sangshen <- df_for_dingxing_analysis_long_sangshen %>%
  select(-c(glycan, peptides, Protein, glycopeptides, patients, `repeat`, value))
df_for_dingxing_analysis_long_sangshen$cancer_glycan=paste(df_for_dingxing_analysis_long_sangshen$cancer_type,df_for_dingxing_analysis_long_sangshen$glycan.type,sep = "-" )
df_for_dingxing_analysis_long_sangshen <- df_for_dingxing_analysis_long_sangshen %>%
  mutate(sites_label = case_when(
    sites == 1 ~ "peptides contain 1 glycan composition",
    sites == 2 ~ "peptides contain 2 glycan composition",
    sites == 3 ~ "peptides contain 3 glycan composition",
    sites == 4 ~ "peptides contain 4 glycan composition",
    sites == 5 ~ "peptides contain 5 glycan composition",
    sites == 6 ~ "peptides contain 6 glycan composition",
    sites > 6 ~ "peptides contain >6 glycan composition",
  ))

df_for_dingxing_analysis_long_sangshen$cancer_glycan_sites=paste(df_for_dingxing_analysis_long_sangshen$cancer_glycan,df_for_dingxing_analysis_long_sangshen$sites_label,sep = "_")

df_for_dingxing_analysis_long_sangshen_count <- df_for_dingxing_analysis_long_sangshen %>%
  dplyr::group_by(cancer_glycan_sites) %>%
  dplyr::summarise(counts = n())
df_for_dingxing_analysis_long_sangshen=merge(df_for_dingxing_analysis_long_sangshen,df_for_dingxing_analysis_long_sangshen_count)
df_for_dingxing_analysis_long_sangshen=df_for_dingxing_analysis_long_sangshen %>% select(-sites)
df_for_dingxing_analysis_long_sangshen=unique(df_for_dingxing_analysis_long_sangshen)
#install.packages("ggalluvial")
library(ggalluvial)
library(ggplot2)
ggplot(df_for_dingxing_analysis_long_sangshen,
       aes(axis1 = cancer_type, axis2 =glycan.type , axis3 = sites_label,
           y = counts)) +
  geom_alluvium(aes(fill = cancer_type)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Glycan Type", "Cancer Type", "Sites")) +
  theme_minimal() +
  labs(title = "Sankey Diagram of Glycan Types, Cancer Types, and Sites",
       x = "Categories", y = "Counts")


proteins_for_go_analysis=unique(df_for_dingxing_analysis_long$Protein)

# 加载必要的库
library(dplyr)
library(VIM)  # 用于KNN填补缺失值
#install.packages("VIM")
# 将NA替换为数值NA，以便后续计算
df_wide[df_wide == "NA"] <- NA
# 计算每个基因在所有样品中的有效值比例，并合并回原始数据框
df_wide <- df_wide %>%
  rowwise() %>%
  mutate(valid_ratio = mean(as.numeric(!is.na(c_across(-contains("Glycopeptides"))))))
class(df_wide$valid_ratio[1])
summary(df_wide$valid_ratio)
# 保留有效值大于50%的基因
df_filtered <- df_wide %>%
  filter(valid_ratio > 0.5) %>%
  dplyr::select(-valid_ratio)
write.csv(df_filtered,"50%有效值保留_knn缺失值填补之前.csv")
# 使用KNN填补缺失值
df_filled <- kNN(df_filtered, variable = setdiff(colnames(df_filtered), "Glycopeptides"))
colnames(df_filled)
df_filled=dplyr::select(df_filled,-ends_with("_imp"))
write.csv(df_filled,file = "填补缺失值.csv")
data_for_cor=df_filled
rownames(data_for_cor) <- data_for_cor$Glycopeptides
data_for_cor=data_for_cor %>% select(-Glycopeptides)
data_for_cor <- data_for_cor[, order(names(data_for_cor))]
# 计算相关性矩阵
correlation_matrix <- cor(data_for_cor, use = "pairwise.complete.obs")
write.csv(correlation_matrix,"样品重复相关性.csv")
# 输出相关性矩阵
print(correlation_matrix)
#install.packages("ggcorrplot")
library(ggplot2)
library(ggcorrplot)
ggcorrplot(correlation_matrix, 
           method = "square", 
           type = "lower", 
           lab = FALSE, 
        #   lab_size = 3, 
           title = "Correlation Matrix Heatmap", 
           ggtheme = theme_classic())


# 将宽数据格式转换为长数据格式
df_long <- df_filled %>%
  pivot_longer(
    cols =  -Glycopeptides,  # 指定要转换的列
    names_to = "Sample",      # 新列的名称，用于存储旧列的名称
    values_to = "Value"       # 新列的名称，用于存储旧列的值
  )
#df_long$Value=log2(df_long$Value)#去log值了
# 绘制抖动散点图加箱线图
p=ggplot(df_long, aes(x = Sample, y = Value)) +
  geom_boxplot() +  # 绘制箱线图，隐藏异常值
  # geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +  # 绘制抖动散点图
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # 调整x轴标签角度
  labs(title = "组间归一化", x = "Sample", y = "Value")+theme_classic()
#ggsave(p,file="组间归一化.pdf",width=25,height=8)
# 将数据框转换为长格式
df_long <- df_filled %>%
  pivot_longer(
    cols = -Glycopeptides,
    names_to = c("Sample", "Repeat"),
    names_pattern = "(.*)_(repeat[0-9]+)"
  )

# 计算每个样品的每个基因的均值
df_mean <- df_long %>%
  group_by(Glycopeptides, Sample) %>%
  summarize(Mean_Expression = mean(value, na.rm = TRUE), .groups = 'drop')
#df_mean$Mean_Expression=log2(df_mean$Mean_Expression)
ggplot(df_mean, aes(x = Sample, y = Mean_Expression)) +
  geom_boxplot() +  # 绘制箱线图，隐藏异常值
 theme_classic() +# geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +  # 绘制抖动散点图
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # 调整x轴标签角度
  labs(title = "Before normalization", x = "Sample", y = "Value")
# 将结果转换为宽格式
df_wide_mean <- df_mean %>%
  pivot_wider(
    names_from = Sample,
    values_from = Mean_Expression
  )

# 查看结果
print(df_wide_mean)
write.csv(df_wide_mean,file = "填补缺失值后.csv")
df_zscore <- df_wide_mean %>%
mutate(across(where(is.numeric) & -Glycopeptides, ~ as.numeric(scale(.))))

df_zscore$Glycopeptides=df_wide_mean$Glycopeptides
df_zscore_long <- df_zscore %>%
  pivot_longer(
    cols = -Glycopeptides,
    names_to = "Sample",      # 新列的名称，用于存储旧列的名称
    values_to = "Value"       # 新列的名称，用于存储旧列的值
  )

ggplot(df_zscore_long, aes(x = Sample, y = Value)) +
  geom_boxplot() +  # 绘制箱线图，隐藏异常值
  theme_classic()+
  # geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +  # 绘制抖动散点图
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  # 调整x轴标签角度为竖直
  labs(title = "After normalization", x = "Sample", y = "Value")  
df_zscore=df_zscore_long%>%pivot_wider(names_from = Sample, values_from = Value)
write.csv(df_zscore,file = "填补缺失值_zscore后.csv")
df=df_zscore

# 加载必要的R包
library(ggplot2)
library(dplyr)
library(tidyr)
#install.packages("ggfortify")
library(ggfortify)
# 提取Glycopeptides列（替代原代码中的第一列）
gene_names <- df$Glycopeptides  # 或 df[["Glycopeptides"]]

# 剩余列名按.分割，并提取注释信息
colnames(df)[-which(names(df) == "Glycopeptides")] <- make.names(colnames(df)[-which(names(df) == "Glycopeptides")])
annotations <- sapply(strsplit(colnames(df)[-which(names(df) == "Glycopeptides")], "\\."), function(x) x[length(x)])
names(df)[-which(names(df) == "Glycopeptides")] <- sapply(strsplit(colnames(df)[-which(names(df) == "Glycopeptides")], "\\."), function(x) x[1])

# 将注释信息添加到数据框
annotation_df <- data.frame(
  Column = colnames(df)[-which(names(df) == "Glycopeptides")],
  Annotation = annotations
)

# 去除Glycopeptides列，保留数值数据
df_numeric <- df %>% select(-Glycopeptides)
# 标准化处理
df_scaled <- scale(df_numeric)
# 执行PCA分析
pca_result <- prcomp(df_scaled, center = TRUE, scale. = TRUE)
# 查看PCA结果
summary(pca_result)
# 将PCA结果转换为数据框
pca_data <- as.data.frame(pca_result$x)
# 添加样本标签
pca_data$Sample <- rownames(pca_data)
# 创建注释信息向量
annotation_vector <- rep(annotations, each = nrow(pca_data))
# 将PC1和PC2的载荷与注释信息结合起来
loadings <- as.data.frame(pca_result$rotation)
loadings$Annotation <- annotations
# 绘制PC1和PC2载荷的注释信息

p2 <- ggplot(loadings, aes(x = PC1, y = PC2, color = Annotation)) +
  geom_point(size = 2) +
  labs(title = "PCA分析", x = "主成分1载荷 (PC1)", y = "主成分2载荷 (PC2)") +
  theme_classic()
print(p2)






