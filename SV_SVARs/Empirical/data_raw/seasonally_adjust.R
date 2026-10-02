setwd("C:/Users/leono/Desktop/Msc Econometrics/thesis/Code/data_raw")

if (!require("seasonal")) install.packages("seasonal")
library(seasonal)


# Read CSV (auto-detects Date column)
data <- read.csv("oecdcrude_unadjusted.csv")

# Convert to time series (monthly)
y <- ts(data$oecdcrude, frequency = 12, start = c(1980, 1))

# Run X-13ARIMA-SEATS
m <- seas(y)

# Extract seasonally adjusted series
data$oecdcrude_sa <- final(m)

plot(final(m))



# Save to CSV
write.csv(data, "oecdcrude_adjusted.csv", row.names = FALSE)