# Glyco FFPE omics 糖肽分析系统 - 依赖包安装脚本
# 版本: v2.0
# 最后更新: 2024年

cat("=== Glyco FFPE omics 依赖包安装脚本 ===\n")
cat("正在检查并安装所需的R包...\n\n")

# 检查R版本
r_version <- R.version.string
cat("当前R版本:", r_version, "\n")

if (as.numeric(R.version$major) < 4 || (as.numeric(R.version$major) == 4 && as.numeric(R.version$minor) < 2)) {
  warning("建议使用R 4.2.0或更高版本以获得最佳性能")
}

# 设置CRAN镜像（可选，提高下载速度）
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# 定义必需的R包列表
required_packages <- c(
  # 核心Shiny框架
  "shiny",
  "shinydashboard", 
  "shinyjs",
  
  # 数据处理和操作
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  
  # 数据可视化
  "ggplot2",
  "ggcorrplot",
  "ggfortify",
  "ggalluvial",
  "ggvenn",
  "ggVennDiagram",
  "ggrepel",
  "ggpubr",
  "patchwork",
  
  # 数据表格和交互
  "DT",
  "plotly",
  "htmlwidgets",
  
  # 热图和相关性分析
  "pheatmap",
  "corrplot",
  
  # 缺失值处理
  "VIM",
  "impute",
  
  # 机器学习和分析
  "factoextra",
  "MASS",
  "e1071",
  
  # ROC分析和模型评估
  "pROC",
  "caret"
)

# Bioconductor包
bioc_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db"
)

# 可选包（用于增强功能）
optional_packages <- c(
  "RColorBrewer",
  "scales",
  "gridExtra",
  "knitr",
  "rmarkdown"
)

# 安装函数
install_if_missing <- function(packages, source = "CRAN") {
  for (package in packages) {
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
      cat("正在安装", package, "...\n")
      tryCatch({
        if (source == "CRAN") {
          install.packages(package, dependencies = TRUE)
        } else if (source == "Bioc") {
          if (!require("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
          }
          BiocManager::install(package, update = FALSE)
        }
        cat("✓", package, "安装成功\n")
      }, error = function(e) {
        cat("✗", package, "安装失败:", e$message, "\n")
      })
    } else {
      cat("✓", package, "已安装\n")
    }
  }
}

# 检查安装状态
check_installation <- function(packages) {
  cat("\n=== 检查安装状态 ===\n")
  failed_packages <- c()
  
  for (package in packages) {
    if (require(package, character.only = TRUE, quietly = TRUE)) {
      cat("✓", package, "已正确安装\n")
    } else {
      cat("✗", package, "安装失败或无法加载\n")
      failed_packages <- c(failed_packages, package)
    }
  }
  
  if (length(failed_packages) > 0) {
    cat("\n警告: 以下包安装失败:\n")
    cat(paste(failed_packages, collapse = ", "), "\n")
    cat("请手动安装这些包或检查网络连接\n")
  } else {
    cat("\n✓ 所有必需包已成功安装!\n")
  }
  
  return(failed_packages)
}

# 主安装流程
cat("=== 开始安装CRAN包 ===\n")
install_if_missing(required_packages, "CRAN")

cat("\n=== 开始安装Bioconductor包 ===\n")
# 检查并安装BiocManager
if (!require("BiocManager", quietly = TRUE)) {
  cat("正在安装BiocManager...\n")
  install.packages("BiocManager")
}
install_if_missing(bioc_packages, "Bioc")

cat("\n=== 开始安装可选包 ===\n")
install_if_missing(optional_packages, "CRAN")

# 检查安装状态
failed_packages <- check_installation(c(required_packages, bioc_packages))

# 加载所有包
cat("\n=== 加载所有包 ===\n")
for (package in c(required_packages, bioc_packages)) {
  tryCatch({
    library(package, character.only = TRUE)
    cat("✓", package, "加载成功\n")
  }, error = function(e) {
    cat("✗", package, "加载失败:", e$message, "\n")
  })
}

# 系统信息
cat("\n=== 系统信息 ===\n")
cat("R版本:", R.version.string, "\n")
cat("操作系统:", Sys.info()["sysname"], Sys.info()["release"], "\n")
cat("内存使用:", round(memory.size() / 1024^2, 2), "MB\n")

# 安装完成提示
if (length(failed_packages) == 0) {
  cat("\n🎉 安装完成! 所有依赖包已成功安装。\n")
  cat("现在可以运行以下命令启动应用:\n")
  cat("source('appchange.R')\n")
} else {
  cat("\n⚠️  安装完成，但部分包安装失败。\n")
  cat("请手动安装失败的包后重新运行此脚本。\n")
}

cat("\n=== 安装脚本结束 ===\n") 