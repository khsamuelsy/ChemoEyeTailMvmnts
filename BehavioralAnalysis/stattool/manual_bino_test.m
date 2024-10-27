% Sample data
% num_successes = 22; % Number of successes observed
% num_trials = 53; % Total number of trials
num_successes = 33; % Number of successes observed
num_trials = 62; % Total number of trials
p0 = 0.5; % Hypothesized proportion of success

% Calculate the p-value for the binomial test
format long
[binocdf(num_successes, num_trials, p0) , 1 - binocdf(num_successes - 1, num_trials, p0)]
p_value = 2 * min(binocdf(num_successes, num_trials, p0), ...
                  1 - binocdf(num_successes - 1, num_trials, p0));

% Determine if we reject the null hypothesis at the 0.05 significance level
alpha = 0.05;
h = p_value < alpha;

% Display the results
disp(['Number of successes: ', num2str(num_successes)]);
disp(['Number of trials: ', num2str(num_trials)]);
disp(['Hypothesized proportion: ', num2str(p0)]);
disp(['Test decision (0 = fail to reject, 1 = reject): ', num2str(h)]);
disp(['p-value: ', num2str(p_value)]);

