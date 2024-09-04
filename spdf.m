% Given CDF data (for demonstration purposes, replace with your actual data)
power_levels = linspace(0, 2, 10000); % example power levels
cdf_values = power_levels / max(power_levels); % example CDF, replace with your data

% Define threshold intervals
thresholds = linspace(0, 2, 10000); % example thresholds, adjust as needed

% Compute the probability of detection for each threshold
prob_detection = 1 - interp1(power_levels, cdf_values, thresholds);

% Form the new CDF
new_cdf = cumsum(prob_detection) / sum(prob_detection);

% Plot the new CDF
figure;
plot(thresholds, new_cdf, 'LineWidth', 1.5);
xlabel('Threshold Level');
ylabel('Cumulative Probability of Detection');
title('CDF of Probability of Detection for Given Thresholds');
grid on;
