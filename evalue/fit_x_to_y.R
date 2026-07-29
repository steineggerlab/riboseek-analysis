library(ggplot2)
library(dplyr)

# Read in the first argv
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Read the data (csv)
data <- read.csv(input_file, header = FALSE)
# Remove the rows if data$V2 == 0
data <- data[data$V2 > 0, ]
data <- data[data$V2 < 1000, ]
# Create new columns for square, exponential, and power of 10
data$square <- data$V2^2
data$exponential <- exp(data$V2)
# Inf to Nan
data$exponential[is.infinite(data$exponential)] <- NA
data$pow10 <- 10^(data$V2)
data$pow10[is.infinite(data$pow10)] <- NA
data$log10 <- log10(data$V2)
data$log10y <- log10(data$V1)

# Fit the second column to the first column and square using a linear model
print("First order")
model <- lm(V1 ~ V2, data = data)
summary(model)
slope <- coef(model)[2]
intercept <- coef(model)[1]
r_squared <- summary(model)$r.squared
cat(paste("Slope:", slope, "\n"))
cat(paste("Intercept:", intercept, "\n"))
cat(paste("R-squared:", r_squared, "\n"))
cat(paste("Equation: y =", slope, "* x +", intercept, "\n"))
data$predicted_first <- predict(model, newdata = data)

# Second order
print("Second order")
model_second <- lm(V1 ~ square + V2, data = data)
summary(model_second)
slope1 <- coef(model_second)[2]
slope2 <- coef(model_second)[3]
intercept <- coef(model_second)[1]
r_squared <- summary(model_second)$r.squared
cat(paste("Slope for x^2:", slope1, "\n"))
cat(paste("Slope for x:", slope2, "\n"))
cat(paste("Intercept:", intercept, "\n"))
cat(paste("R-squared:", r_squared, "\n"))
cat(paste("Equation: y =", slope1, "* x^2 +", slope2, "* x +", intercept, "\n"))
data$predicted_second <- predict(model_second, newdata = data)

# Exponential
print("Exponential")
model_exp <- lm(V1 ~ exponential, data = data, na.action = na.exclude)
summary(model_exp)
slope <- coef(model_exp)[2]
intercept <- coef(model_exp)[1]
r_squared <- summary(model_exp)$r.squared
cat(paste("Slope:", slope, "\n"))
cat(paste("Intercept:", intercept, "\n"))
cat(paste("R-squared:", r_squared, "\n"))
cat(paste("Equation: y =", slope, "* exp(x) +", intercept, "\n"))
data$predicted_exp <- predict(model_exp, newdata = data)

# Power of 10
print("Power of 10")
model_10 <- lm(V1 ~ pow10, data = data, na.action = na.exclude)
summary(model_10)
slope <- coef(model_10)[2]
intercept <- coef(model_10)[1]
r_squared <- summary(model_10)$r.squared
cat(paste("Slope:", slope, "\n"))
cat(paste("Intercept:", intercept, "\n"))
cat(paste("R-squared:", r_squared, "\n"))
cat(paste("Equation: y =", slope, "* 10^x +", intercept, "\n"))
data$predicted_10 <- predict(model_10, newdata = data)

# Log10
print("Log10")
model_log10 <- lm(log10y ~ log10, data = data, na.action = na.exclude)
summary(model_log10)
slope <- coef(model_log10)[2]
intercept <- coef(model_log10)[1]
r_squared <- summary(model_log10)$r.squared
cat(paste("Slope:", slope, "\n"))
cat(paste("Intercept:", intercept, "\n"))
cat(paste("R-squared:", r_squared, "\n"))
cat(paste("Equation: y =", slope, "* log10(x) +", intercept, "\n"))
data$predicted_log10 <- 10^predict(model_log10, newdata = data)

print(head(data))

# Make a plot of the data and the fitted lines
# plot X axis and Y axis as log scale
# x and y ranging to 1000
# Also draw y=x
df <- data.frame(
  value = 1:1000
)

# V1 is reported evalue, V2 is observed evalue
p <- ggplot(data, aes(x = V1, y = V2)) +
  geom_point(color = 'black') +
  geom_line(aes(y = predicted_first), color = 'blue', size = 1, linetype = "dashed") +
  geom_line(aes(y = predicted_second), color = 'red', size = 1, linetype = "dashed") +
  geom_line(aes(y = predicted_exp), color = 'green', size = 1, linetype = "dashed") +
  geom_line(aes(y = predicted_10), color = 'purple', size = 1, linetype = "dashed") +
  geom_line(aes(y=predicted_log10), color = 'orange', size = 1, linetype = "dashed") +
  labs(title = 'Fitting y to x with Different Models',
       x = 'X values',
       y = 'Y values') +
  theme_minimal() +
  scale_y_continuous(limits = c(min(data$V2, na.rm = TRUE) * 0.9, 1000)) +
  scale_y_log10() +
  scale_x_continuous(limits = c(min(data$V1) * 0.9, 1000)) +
  scale_x_log10() +
  theme(legend.position = "top") +
  guides(color = guide_legend(title = "Models")) +
  annotate("text", x = Inf, y = Inf, label = paste("First order R²:", round(summary(model)$r.squared, 4)), hjust = 1.1, vjust = 2, color = "blue", size = 3) +
  annotate("text", x = Inf, y = Inf, label = paste("Second order R²:", round(summary(model_second)$r.squared, 4)), hjust = 1.1, vjust = 3.5, color = "red", size = 3) +
  annotate("text", x = Inf, y = Inf, label = paste("Exponential R²:", round(summary(model_exp)$r.squared, 4)), hjust = 1.1, vjust = 5, color = "green", size = 3) +
  annotate("text", x = Inf, y = Inf, label = paste("Power of 10 R²:", round(summary(model_10)$r.squared, 4)), hjust = 1.1, vjust = 6.5, color = "purple", size = 3) +
  annotate("text", x = Inf, y = Inf, label = paste("Log10 R²:", round(summary(model_log10)$r.squared, 4)), hjust = 1.1, vjust = 8, color = "orange", size = 3) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")

# Save the plot
ggsave(output_file, plot = p, width = 8, height = 6)
print(paste("Plot saved to", output_file))