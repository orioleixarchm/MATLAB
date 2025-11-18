%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Advanced Econometrics                                                   %
% Sample Bias and Bootstrap Correction in Finite Samples for AR(1) Models %
% Oriol Eixarch Mejías r0872954                                           %
% 21/09/2025                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cleanup
clear
close all
clc

% Parameters
rng(111)                            % Setting seed
rho_vec = [0.3 0.7];                % True rho values (one low, one high)
T_vec = [5 10 15 20 30 45 60 75];   % Sample sizes      
R = 2000;                           % Monte Carlo replications
B = 1000;                           % Bootstrap replications
buff = 100;                         % Burn-in observations to remove initial effect first value (the initial contribution dies out at a rate rho)
alpha_ci = 0.05;                    % 95% Confidence interval

% Precreation of Figures (A = simple CI, B = bootstrap-based CI)
fig1 = figure('Name','Distributions: rho-hat vs rho^{BC} (Option A)');
fig2 = figure('Name','SE vs SD comparison (Option A)');
fig3 = figure('Name','Distributions: rho-hat vs rho^{BC} (Option B)');
fig4 = figure('Name','SE vs SD comparison (Option B)');

tl1  = tiledlayout(fig1,2,4,'Padding','compact','TileSpacing','compact');
tl2  = tiledlayout(fig2,2,4,'Padding','compact','TileSpacing','compact');
tl3  = tiledlayout(fig3,2,4,'Padding','compact','TileSpacing','compact');
tl4  = tiledlayout(fig4,2,4,'Padding','compact','TileSpacing','compact');

sgtitle(tl1,sprintf('Distributions of $\\hat{\\rho}$ and $\\hat{\\rho}^{BC}$ (Option A, $\\rho = %.1f$)', rho_vec(1)),'Interpreter','latex')
sgtitle(tl2,sprintf('Mean OLS SE vs Empirical SD vs Mean Bootstrap SE (Option A, $\\rho = %.1f$)',rho_vec(1)),'Interpreter','latex')
sgtitle(tl3,sprintf('Distributions of $\\hat{\\rho}$ and $\\hat{\\rho}^{BC}$ (Option B, $\\rho = %.1f$)', rho_vec(2)),'Interpreter','latex')
sgtitle(tl4,sprintf('Mean OLS SE vs Empirical SD vs Mean Bootstrap SE (Option B, $\\rho = %.1f$)',rho_vec(2)),'Interpreter','latex')

% Rho Loop
for rho = rho_vec
    p = 0; % plot index
    for T = T_vec
        p = p + 1; %update plot index

        % Preallocate
        rho_hat_all = zeros(R,1);
        rho_bc_all  = zeros(R,1);
        ols_se_all  = zeros(R,1);
        boot_se_all = zeros(R,1);
        bc_se_all   = zeros(R,1);          % Used only for Option B confidence interval

        co_ols = zeros(R,1);
        co_bc_A = zeros(R,1);
        co_bc_B = zeros(R,1);

        ci_ols_store = zeros(R,2);
        ci_bc_store_A = zeros(R,2);
        ci_bc_store_B = zeros(R,2);

        % Constants
        k  = 1;                                
        df = T - k - 1;   %residual degrees of freedom T-nºparameters - 1 unusable pair (first pair)                  
        crit = tinv(1 - alpha_ci/2, df);   %tinv(p, df) inverse cdf of the Student's t distribution / the p-th quantile.      

        % Monte Carlo Loop
        for r=1:R
            % Simulate AR(1) process
            y_ph = zeros(T+buff,1);
            u_ph = randn(T+buff,1);
            y_ph(1) = 2;
            for i = 2:(T+buff)
                y_ph(i) = rho*y_ph(i-1) + u_ph(i);
            end
            y_ph = y_ph(end-T+1:end);
            
            y = y_ph(2:end);
            x = y_ph(1:end-1);

            % OLS estimation 
            rho_hat = ((x'*x)^(-1))*(x'*y);
            u = y - rho_hat*x;
            u_c = u - mean(u);   %bootstrapping requires errors to have mean zero               
            sigma2_hat = (u'*u) / (T-2); %For ols we require raw errors even if they are not mean 0 (In finite samples) otherwise we underestimate the variance u'u>u_c'U_c       
            ols_se = sqrt(sigma2_hat / (x'*x)); 
            
            % CI (OLS)
            ci_ols = [rho_hat - crit*ols_se, rho_hat + crit*ols_se];
            ci_ols_store(r,:) = ci_ols;
            co_ols(r) = (rho >= ci_ols(1)) && (rho <= ci_ols(2));  %1 if true rho is inside CI 0 otherwises (OLS)

            % Bootstrap (Residual)
            rho_star = zeros(B,1);
            for i = 1:B
                u_star_ph = u_c(randi(T-1,T-1,1));  
                y_star_ph = zeros(T,1);
                y_star_ph(1) = y(1);
                for t = 2:T
                    y_star_ph(t) = rho_hat*y_star_ph(t-1) + u_star_ph(t-1);
                end
                y_star = y_star_ph(2:end);
                x_star = y_star_ph(1:end-1);
                rho_star(i) = ((x_star'*x_star)^(-1))*(x_star'*y_star);
            end

            % Bias Correction 
            mean_rho = mean(rho_star);
            boot_se_s = std(rho_star,0);
            bias = mean_rho - rho_hat;
            rho_bc = rho_hat - bias;

            % Option A: use OLS SE 
            se_A = ols_se;
            ciA = [rho_bc - crit*se_A, rho_bc + crit*se_A];

            % Option B: use bootstrap SE of BC estimator 
            rho_bc_star = 2*rho_star - mean_rho; %Equivalent to Boostrap Replication%
            se_B = std(rho_bc_star,0);
            ciB = [rho_bc - crit*se_B, rho_bc + crit*se_B];

            % Store results
            ci_bc_store_A(r,:) = ciA;
            ci_bc_store_B(r,:) = ciB;
            co_bc_A(r) = (rho >= ciA(1)) && (rho <= ciA(2));
            co_bc_B(r) = (rho >= ciB(1)) && (rho <= ciB(2));

            rho_hat_all(r) = rho_hat;
            rho_bc_all(r)  = rho_bc;
            ols_se_all(r)  = ols_se;
            boot_se_all(r) = boot_se_s;
            bc_se_all(r)   = se_B;
        end

        % Aggregate results
        mean_rho_hat = mean(rho_hat_all); %Average OLS estimator across Monte Carlo estimations
        mean_rho_bc  = mean(rho_bc_all); %Average bias-corrected estimator across Montecarlo estimations
        bias_rho_hat = mean_rho_hat - rho; %Montecarlo Bias of OLS
        bias_rho_bc  = mean_rho_bc - rho; %Montecarlo Bias of the bias-corrected estimator
        se_rho_hat   = std(rho_hat_all);  %Standard deviation of the OLS
        mean_se_ols  = mean(ols_se_all); %Mean SE OLS
        mean_se_boot_all = mean(boot_se_all); %Mean bootstrap SE across Montecarlo replications
        mean_se_bc_B = mean(bc_se_all);
        diff = mean_rho_bc - mean_rho_hat;

        % Coverage rates
        %Percentage of times estimator within the CI
        coverage_rate_ols  = mean(co_ols); 
        coverage_rate_bc_A = mean(co_bc_A);
        coverage_rate_bc_B = mean(co_bc_B);

        % Summary Table
        summ = table(rho, mean_rho_hat, mean_rho_bc, bias_rho_hat, bias_rho_bc, diff, ...
                     coverage_rate_ols, coverage_rate_bc_A, coverage_rate_bc_B, ...
                     mean(ci_ols_store,1), mean(ci_bc_store_A,1), mean(ci_bc_store_B,1), ...
                     'VariableNames',{'TrueRho','MeanOLS','MeanBC','BiasOLS','BiasBC','BCvsOLS_Diff', ...
                     'Coverage_OLS_95','Coverage_BC_A_95','Coverage_BC_B_95', ...
                     'Mean_CI_OLS','Mean_CI_BC_A','Mean_CI_BC_B'});

        fprintf('\n================ Summary for T = %d and rho = %.1f ================\n', T, rho);
        disp(summ);

        % Plotting
        if rho == rho_vec(1)
            nexttile(tl1,p)
            histogram(rho_hat_all,'Normalization','pdf','FaceAlpha',0.35,'EdgeColor','none'); hold on
            histogram(rho_bc_all,'Normalization','pdf','FaceAlpha',0.35,'EdgeColor','none');
            xline(rho,'-','LineWidth',1.5);
            xline(mean_rho_hat,':','LineWidth',1.5);
            xline(mean_rho_bc,'-.','LineWidth',1.5);
            title(sprintf('T = %d',T))
            xlabel('\rho estimate'); ylabel('Density');
            legend({'$\hat{\rho}$ (OLS)', '$\hat{\rho}^{BC}$', 'true $\rho$', ...
                    'mean $\hat{\rho}$', 'mean $\hat{\rho}^{BC}$'}, ...
                    'Interpreter','latex','Location','best');
            grid on

            nexttile(tl2,p)
            bar_vals = [mean_se_ols, se_rho_hat, mean_se_boot_all,mean_se_bc_B];
            bar(categorical({'Scenario'}), bar_vals);
            ylabel('Standard Errors and SD');
            title(sprintf('T = %d',T))
            legend({'Mean OLS SE','Empirical SD($\hat{\rho}$)','Mean bootstrap SE','Mean SD($\hat{\rho}^{BC}$)'}, ...
                    'Interpreter','latex','Location','northwest');
            grid on
        else
            nexttile(tl3,p)
            histogram(rho_hat_all,'Normalization','pdf','FaceAlpha',0.35,'EdgeColor','none', ...
                'FaceColor',[0.0 0.2 0.4]); hold on
            histogram(rho_bc_all,'Normalization','pdf','FaceAlpha',0.35,'EdgeColor','none', ...
                'FaceColor',[0.6000 0.2000 0.2000]);
            xline(rho,'-','LineWidth',1.5);
            xline(mean_rho_hat,':','LineWidth',1.5);
            xline(mean_rho_bc,'-.','LineWidth',1.5);
            title(sprintf('T = %d',T))
            xlabel('\rho estimate'); ylabel('Density');
            legend({'$\hat{\rho}$ (OLS)', '$\hat{\rho}^{BC}$', 'true $\rho$', ...
                    'mean $\hat{\rho}$', 'mean $\hat{\rho}^{BC}$'}, ...
                    'Interpreter','latex','Location','best');
            grid on

            nexttile(tl4,p)
            bar_vals = [mean_se_ols, se_rho_hat, mean_se_boot_all,mean_se_bc_B];
            bb = bar(categorical({'Scenario'}), bar_vals);
            ylabel('Standard Errors and SD');
            bb(1).FaceColor = [0.0 0.2 0.4];
            bb(2).FaceColor = [0.6000 0.2000 0.2000];
            bb(3).FaceColor = [0.8500 0.6000 0.2000];
            title(sprintf('T = %d',T))
            legend({'Mean OLS SE','Empirical SD($\hat{\rho}$)','Mean bootstrap SE', 'Mean SD($\hat{\rho}^{BC}$)'}, ...
                    'Interpreter','latex','Location','northwest');
            grid on
        end
    end
end

%% OPTION A / B EXPLANATION
% Option A :
%   Explanation: We consider that the bias correction is a deterministic
%   shift of the center and thus we can still use rho's varaince (ols). It is
%   simpler but comes at the cost of ignoring the extra uncertainty from
%   estimating the bias.
%   create a bootstrap analog for our BC estimation, we compute its se.

%   Why is our coverage lower for BC than for OLS?
%   OLS se is lower than the actual BC se, which means that our CI are
%   naively narrow wich is why the coverage is also lower (more points are
%   outside the bounds). 
%   Treats bias correction as a deterministic shift.
%   Uses OLS SE for CI, ignores uncertainty in bias estimate.
%   Pros: simple; Cons: coverage may be too low (CIs too narrow).

% Option B (Bootstrap-based):
%   Explanation: We manually try to use the formula found in Efron, 1979  to
%   create a bootstrap analog for our BC estimation, we compute its se.
%   Closer to reality, se becomes very big due to the bias inclusion which is 
%   essentially a random variable,thus CI becames huge and thus more concervative
%   (reaching a coverage of 100% sometimes); on top of that, it adds complexity 
%   to the code and interpretation.

%   Why are our BC se significantly higher than for OLS (why are our BC CI wider than those of OLS)?
%   rho^{BC} is computed by substracting the bias=ρ^​∗​−ρ^ to rho_hat, given that
%   the bias is a random variable (it depends on bootstrap draws), adds extra variability. 
%   We face a clear trade-off increased variability in exchange of a less biase estimator, 
%   that is why CI for the BC estimator are significantly widder.

%   Bias of rho when regressing y_t-1 on y_ t => E[ρ_hat−ρ]≈−(1+3ρ)/T, larger values of rho (sample ie unchanged) lead to larger bias, as the sample grows, the bias fades. 
%   Formula extracted from existing literature e.g: Kendall 1954, Nickell 1958
%   Accounts for bias estimation variability via bootstrap SEs of the BC estimator.
%   rho_bc_star = 2*rho_star - mean(rho_star).
%   Pros: more accurate uncertainty; Cons: CIs wider, sometimes overly conservative.
