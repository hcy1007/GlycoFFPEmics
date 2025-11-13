# 测试组合图表功能
library(ggplot2)
library(patchwork)

# 创建测试数据
set.seed(123)
test_data <- data.frame(
  x = rnorm(100),
  y = rnorm(100),
  group = sample(c("A", "B", "C"), 100, replace = TRUE)
)

# 创建四个测试图表
plot1 <- ggplot(test_data, aes(x = x, y = y, color = group)) +
  geom_point() +
  ggtitle("Plot 1") +
  theme_minimal()

plot2 <- ggplot(test_data, aes(x = group, y = x, fill = group)) +
  geom_boxplot() +
  ggtitle("Plot 2") +
  theme_minimal()

plot3 <- ggplot(test_data, aes(x = x)) +
  geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7) +
  ggtitle("Plot 3") +
  theme_minimal()

plot4 <- ggplot(test_data, aes(x = y)) +
  geom_density(fill = "orange", alpha = 0.7) +
  ggtitle("Plot 4") +
  theme_minimal()

# 创建2x2组合图表
combined_plot <- (plot1 + plot2) / (plot3 + plot4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Test Combined Plot",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

# 显示组合图表
print(combined_plot)

# 保存组合图表
ggsave("test_combined_plot.png", combined_plot, width = 12, height = 10, dpi = 300)

cat("测试完成！组合图表已保存为 test_combined_plot.png\n") 