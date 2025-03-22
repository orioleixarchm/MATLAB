% Filename: sticky_prices_wages_loop_rules_B_matlab.m

% List of iterative variables
Phi_p_1  = 0:0.1:1.4;
Phi_w_1  = 0:0.1:1.4; 
Phi_p_2  = 1.5 * ones(1,15);
Phi_w_2  = 1.5 * ones(1,15); 
Phi_p    = [Phi_p_1, Phi_p_2];
Phi_w    = [Phi_w_2, Phi_w_1];
N        = length(Phi_p);
save pars1 Phi_p Phi_w N

% Loss values
L = zeros(N,1);

% Looping and solving the model for each parameter value
for j=1:N
    load pars1
    curr_phi_p = Phi_p(j);
    curr_phi_w = Phi_w(j);
    save pars2 curr_phi_p curr_phi_w j
    dynare sticky_prices_wages_loop_rules_B;
    load pars2
    % Time series extraction
    pi_p_sim = oo_.endo_simul(strmatch('pi_p', M_.endo_names, 'exact'), :);
    pi_w_sim = oo_.endo_simul(strmatch('pi_w', M_.endo_names, 'exact'), :);
    y_tilde_sim = oo_.endo_simul(strmatch('y_tilde', M_.endo_names, 'exact'), :);
    var_pi_p = var(pi_p_sim);
    var_pi_w = var(pi_w_sim);
    var_y_tilde = var(y_tilde_sim);
    % Parameters extraction
    alpha = M_.params(strmatch('alpha', M_.param_names, 'exact'));
    sigma = M_.params(strmatch('sigma', M_.param_names, 'exact'));
    phi = M_.params(strmatch('phi', M_.param_names, 'exact'));
    epsilon_p = M_.params(strmatch('epsilon_p', M_.param_names, 'exact'));
    epsilon_w = M_.params(strmatch('epsilon_w', M_.param_names, 'exact'));
    theta_p = M_.params(strmatch('theta_p', M_.param_names, 'exact'));
    theta_w = M_.params(strmatch('theta_w', M_.param_names, 'exact'));
    beta = M_.params(strmatch('beta', M_.param_names, 'exact'));
    lambda_p = (1 - theta_p) * (1 - beta * theta_p) / (theta_p * ((1 - alpha) / (1 - alpha + alpha * epsilon_p)));
    lambda_w = (1 - theta_w) * (1 - beta * theta_w) / (theta_w * (1 + epsilon_w * phi));
    % Loss function
    L(j,1) = 0.5 * ((sigma + (phi + alpha) / (1 - alpha)) * var_y_tilde + (epsilon_p / lambda_p) * var_pi_p + (epsilon_w * (1 - alpha) / lambda_w) * var_pi_w);
    save Loss_loop L
end

% Create a table with the calculated moments
for a=1:N
    Values{a} = [sprintf('Phi_p = %.2f, Phi_w = %.2f',Phi_p(a), Phi_w(a))];
end

% Combine into a table
results_table = table(Values, L);

% Display the table
disp('Loss Table:');
disp(results_table);
