% Filename: sticky_prices_wages.mod

% Declaring endogenous variables
var pi_p pi_w y_tilde w_tilde wn_t r_nt i_t a_t z_t nu_t;

% Declaring exogenous shocks
varexo epsilon_a epsilon_z epsilon_nu;

% Declaring model parameters
parameters alpha beta sigma phi epsilon_p epsilon_w theta_p theta_w rho_nu rho_a rho_z sigma_nu sigma_a sigma_z phi_p phi_w;

alpha = 0.25;         % Capital share
beta = 0.99;          % Quarterly discount factor
sigma = 2.0;          % CRRA utility coefficient
phi = 2.0;            % Inverse of Frisch elasticity of labor supply
epsilon_p = 9;        % Elasticity of substitution (goods)
epsilon_w = 4.5;      % Elasticity of substitution (labor)
theta_p = 0.75;       % Price rigidity
theta_w = 0.75;       % Wage rigidity
rho_nu = 0.5;         % Persistence of monetary policy shock
rho_a = 0.9;          % Persistence of technology shock
rho_z = 0.5;          % Persistence of preference shock
sigma_nu = 0.01;      % Standard deviation of monetary policy shock
sigma_a = 0.01;       % Standard deviation of technology shock
sigma_z = 0.02;       % Standard deviation of preference shock
phi_p = 1.5;          % Taylor rule coefficient for price inflation
phi_w = 0.0;          % Taylor rule coefficient for wage inflation

% Complex parameters
parameters psi_wa psi_ya kappa_p lambda_p kappa_w lambda_w;

% Derived parameters
lambda_p = (1 - theta_p) * (1 - beta * theta_p) / (theta_p * ((1 - alpha) / (1 - alpha + alpha * epsilon_p)));
lambda_w = (1 - theta_w) * (1 - beta * theta_w) / (theta_w * (1 + epsilon_w * phi));
kappa_p = alpha * lambda_p / (1 - alpha);
kappa_w = lambda_w * (sigma + phi / (1 - alpha));
psi_wa = (1 - alpha * phi) / (1 - alpha);
psi_ya = (1 + phi) / (sigma * (1 - alpha) + phi + alpha);

% Model equations
model;
    % Price Phillips curve
    pi_p = beta * pi_p(+1) + kappa_p * y_tilde + lambda_p * w_tilde;

    % Wage Phillips curve
    pi_w = beta * pi_w(+1) + kappa_w * y_tilde - lambda_w * w_tilde;

    % Wage gap dynamics
    w_tilde = w_tilde(-1) + pi_w - pi_p - wn_t + wn_t(-1);

    % Dynamic IS curve
    y_tilde = -(1/sigma) * (i_t - pi_p(+1) - r_nt) + y_tilde(+1);

    % Taylor rule
    i_t = phi_p * pi_p + phi_w * pi_w + nu_t;

    % Natural real wage
    wn_t = psi_wa * a_t;

    % Natural real interest rate
    r_nt = -sigma * (1 - rho_a) * psi_ya * a_t + (1 - rho_z) * z_t;

    % Shocks
    a_t = rho_a * a_t(-1) + epsilon_a;
    z_t = rho_z * z_t(-1) + epsilon_z;
    nu_t = rho_nu * nu_t(-1) + epsilon_nu;
end;

% Steady state block
% steady_state_model;
%     a_t = 0;
%     z_t = 0;
%     nu_t = 0;
%     pi_p = 0;
%     pi_w = 0;
%     w_tilde = 0;
%     y_tilde = 0;
%     r_nt = 0;
%     i_t = 0;
%     wn_t = 0;
% end;

% Compute steady state
steady;

% Shocks definition
shocks;
    var epsilon_a = sigma_a^2;
    var epsilon_z = sigma_z^2;
    var epsilon_nu = sigma_nu^2;
end;

% Simulation settings
stoch_simul(periods=10000, order=1, irf=20);
