% Load Dynare results
dynare sticky_prices_wages_I.mod;

% Extract simulated variables
pi_p_sim = oo_.endo_simul(strmatch('pi_p', M_.endo_names, 'exact'), :);
pi_w_sim = oo_.endo_simul(strmatch('pi_w', M_.endo_names, 'exact'), :);
y_tilde_sim = oo_.endo_simul(strmatch('y_tilde', M_.endo_names, 'exact'), :);

% Calculate standard deviations
std_pi_p = std(pi_p_sim);
std_pi_w = std(pi_w_sim);
std_y_tilde = std(y_tilde_sim);

% Create a table with the calculated moments
Variable = [' Output Gap (y_tilde) '; 'Price Inflation (pi_p)'; 'Wage Inflation (pi_w) '];
Standard_Deviation = [std_y_tilde; std_pi_p; std_pi_w];

% Combine into a table
results_table = table(Variable, Standard_Deviation);

% Display the table
disp('Calculated Moments Table:');
disp(results_table);
