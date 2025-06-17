library(ggplot2)

# Sample data
df <- data.frame(
  x = 1:10,
  y1 = rnorm(10, 10, 2),
  y2 = rnorm(10, 100, 20)
)

# Create the plot
ggplot(df, aes(x = x)) +
  geom_line(aes(y = y1), color = "blue") + 
  geom_line(aes(y = y2 / 10), color = "red") +  # Adjust y2 to match the scale of y1
  scale_y_continuous(
    name = "Primary Y Axis (y1)",
    sec.axis = sec_axis(~ . * 10, name = "Secondary Y Axis (y2)")  # Adjust back to original scale
  ) +
  theme_minimal() +
  labs(
    title = "Plot with Two Different Y-Axis Scales",
    x = "X Axis"
  )
