%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Advanced Econometrics                                                   %
% Sample Bias and Bootstrap Correction in Finite Samples for AR(1) Models %
% Kaat Verbist         r0840862                                           %
% Oriol Eixarch Mejías r0872954                                           %
% 21/09/2025                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cleanup
clear
close all
clc

% Parameters
rng(111)                            %Setting seed
rho_vec = [0.3 0.7];                %True parameters of rho; one low one high to ppreciate differences in bias (bias higher for hiher rho)
T_vec = [5 10 15 20 30 45 60 75];   %Sample Sizes      
R = 2000;                           %Montecarlo replications
B = 1000;                           %Bootstrap replications
buff = 100;                         %Buffer Observations to remove the effect of the first value (the initial contribution dies out at a rate rho)
alpha_ci = 0.05;                    %95% Conf interval

co_ols = zeros(R,1);                %Result Grid OLS estimator 1 i within CI 0 otherwise
co_bc = zeros(R,1);                 %Result Grid BC estimator 1 i within CI 0 otherwise

ci_ols_store  = zeros(R,2);         %OLS CI Endpoints
ci_bc_store   = zeros(R,2);         %BC CI Endpoints

% Precreation of Subplots
fig1 = figure('Name','Distributions: rho-hat vs rho^{BC} A');
fig2 = figure('Name','SE vs SD comparison A');
fig3 = figure('Name','Distributions: rho-hat vs rho^{BC} B');
fig4 = figure('Name','SE vs SD comparison B');
tl1  = tiledlayout(fig1,2,4,'Padding','compact','TileSpacing','compact');
tl2  = tiledlayout(fig2,2,4,'Padding','compact','TileSpacing','compact');
tl3  = tiledlayout(fig3,2,4,'Padding','compact','TileSpacing','compact');
tl4  = tiledlayout(fig4,2,4,'Padding','compact','TileSpacing','compact');
sgtitle(tl1,sprintf('Distributions of $\\hat{\\rho}$ and $\\hat{\\rho}^{BC}$ ($\\rho = %.1f$)', rho_vec(1)),'Interpreter','latex')
sgtitle(tl2,sprintf('Mean OLS SE vs Empirical SD vs Mean Bootstrap SE ($\\rho = %.1f$)',rho_vec(1)),'Interpreter','latex')
sgtitle(tl3,sprintf('Distributions of $\\hat{\\rho}$ and $\\hat{\\rho}^{BC}$ ($\\rho = %.1f$)', rho_vec(2)),'Interpreter','latex')
sgtitle(tl4,sprintf('Mean OLS SE vs Empirical SD vs Mean Bootstrap SE ($\\rho = %.1f$)',rho_vec(2)),'Interpreter','latex')

% Rho Loop
for rho = rho_vec
    p = 0;          %Plot index
    for T = T_vec
        p = p+1;    %Updating plot index
    
        % Outer Montecarlo Loop
        rho_hat_all = zeros(R,1);
        rho_bc_all = zeros(R,1);
        ols_se_all = zeros(R,1);
        boot_se_all = zeros(R,1);
        %bc_se_all = zeros(R,1);                %%Used in option B for SE and CI computations%%           
        
        k  = 1;                                 %Only parameter to estimate AR(1)
        df = T - k - 1;                         %residual degrees of freedom = T-nºparameters - 1 unusable pair (first pair)
        crit = tinv(1 - alpha_ci/2, df);        %tinv(p, df) inverse cdf of the Student’s t distribution / the p-th quantile.

        for r=1:R
            % Generation of our AR(1) Process
            y_ph = zeros(T+buff,1);
            u_ph = randn(T+buff,1);
            y_ph(1) = 2;
            for i = 2:(T+buff)
                y_ph(i) = rho*y_ph(i-1) + u_ph(i);
            end
            y_ph = y_ph(end-T+1:end);
            
            % OLS Simulation
            y = y_ph(2:end);
            x = y_ph(1:end-1);
            rho_hat = ((x'*x)^(-1))*(x'*y);
            u = y - rho_hat*x;
            u_c = u - mean(u);                  %For bootstrapping we require errors to have mean 0
            sigma2_hat = (u'*u) / (T-2);        %For ols we require raw errors even if they are not mean 0 (In finite samples) otw we underestimate the variance u'u>u_c'U_c
            ols_se = sqrt(sigma2_hat / (x'*x)); 
            
            ci_ols = [rho_hat - crit*ols_se,  rho_hat + crit*ols_se];   %Construction of CI intervals (OLS)
            ci_ols_store(r,:) = ci_ols;                                 %Store CI endpoints for this replications (OLS)
            co_ols(r) = (rho >= ci_ols(1)) && (rho <= ci_ols(2));       %1 if true rho is inside CI 0 otherwises (OLS)

            % Bootstrapping
            rho_star = zeros(B,1);
            for i = 1:B
                u_star_ph = u_c(randi(T-1,T-1,1));  %residuals from 2:T (Need to use t-1 in the inner loop)
                y_star_ph = zeros(T,1);
                y_star_ph(1) = y(1);
                for t = 2:T
                    y_star_ph(t) = rho_hat*y_star_ph(t-1) + u_star_ph(t-1); %u_c(randi(T-1,T-1,1)) creates a vector of length 𝑇−1, notT, need to align backward”.
                end
                y_star = y_star_ph(2:end);
                x_star = y_star_ph(1:end-1);
                rho_star(i) = ((x_star'*x_star)^(-1))*(x_star'*y_star);
            end
            
            % Saving Rho_hat, Rho_bc, standard error of Rho* (Bootstrapped)
            mean_rho = mean(rho_star);
            boot_se_s = std(rho_star,0);                            %Sample standard error (Bessel correction applied)
            bias = mean_rho - rho_hat;                              %Bias appears since, in small samples, the corr between y_t-1 and u_t does not vanish / is non-zero
            rho_bc = rho_hat - bias;                                %Corrected estimator rho
       
            %rho_bc_star   = 2*rho_star - mean_rho;                 %%Used in option B for SE and CI computations IN Boostrap Replication%%
                                                                    %bootstrap analogs of rho_bc, rho_bc = rho_hat - bias: rho_hat - (mean_rho - rho_hat) = 2*rho_hat - mean_rho
            %se_bc_boot    = std(rho_bc_star, 0);                   %%Used in option B for SE and CI computations%%
                                                                    %Se of our Rho BC estimator
            se_bc_boot = ols_se;
            %bc_se_all(r) = se_bc_boot;                             %%Used in option B for SE and CI computations%%
            
            % Confidence intervals for Rho BC (option A)
            ci_bc = [rho_bc - crit.*se_bc_boot, rho_bc  + crit.*se_bc_boot];
            ci_bc_store(r,:) = ci_bc;
            co_bc= (rho >= ci_bc_store(:,1)) & (rho <= ci_bc_store(:,2));

            rho_hat_all(r) = rho_hat;
            rho_bc_all(r) = rho_bc;
            boot_se_all(r) = boot_se_s;
            ols_se_all(r) = ols_se;                                 
        end
        
        % Confidenec intervals for Rho BC (option B)
        %ci_bc_store = [rho_bc_all - crit.*bc_se_all,  rho_bc_all + crit.*bc_se_all];   %%Used in option B for SE and CI computations%%   
                                                                                        %Construction of CI intervals (BC) + Store CI endpoints for each replication (BC)
        %co_bc= (rho >= ci_bc_store(:,1)) & (rho <= ci_bc_store(:,2));                  %%Used in option B for SE and CI computations%%
                                                                                        %1 if true rho is inside CI 0 otherwises (BC)
        
        % Computing main interest statistics
        mean_rho_hat = mean(rho_hat_all);       %Average OLS estimator across MC estimations E(rho_hat)
        bias_rho_hat = mean_rho_hat - rho;      %Montecarlo Bias of OLS
        se_rho_hat = std(rho_hat_all);          %Standard deviation of the OLS estimates across MC replications

        mean_rho_bc = mean(rho_bc_all);         %Average bias-corrected estimator across MC estimations E(rho_bc)
        bias_rho_bc = mean_rho_bc - rho;        %Montecarlo Bias of the bias-corrected estimator
        %se_rho_bc = std(rho_bc_all);           %%Used in option B for SE and CI computations%%
                                                %Standard deviation of the BC estimates across MC replications

        mean_se_boot_all = mean(boot_se_all);   %Average bootstrap standard error across MC replications.
        mean_se_ols = mean(ols_se_all);         %Mean Standard Error of the OLS estimator.
        %mean_se_bc = mean(bc_se_all);          %%Used in option B for SE and CI computations%%
                                                %Mean Standard Error of the BC estimator.

        diff = mean_rho_bc-mean_rho_hat;
        
        coverage_rate_ols = mean(co_ols);       %Percentage of times OLS estimator within the CI
        coverage_rate_bc  = mean(co_bc);        %Percentage of times BC estimator within the CI

        % Display Table construction
        summ = table(rho, mean_rho_hat, mean_rho_bc, bias_rho_hat, ...
                     bias_rho_bc, diff, ...
                     coverage_rate_ols, coverage_rate_bc, ...
                     mean(ci_ols_store,1), mean(ci_bc_store,1), ...
                     'VariableNames',{'True Rho', 'Mean OLS Rho', ...
                     'Mean Bias-Corrected Rho','Bias OLS Rho',  ...
                     'Bias MC Bias-Corrected Rho','Bias-Corrected Rho - Mean OLS Rho', ...
                     'Coverage OLS (95%)','Coverage BC (95%)','Mean of CI Endpoints OLS', ...
                     'Mean of CI Endpoints BC'});
    
        %%Auxiliar Var used in option B%% se_rho_hat, mean_se_ols, mean_se_boot_all, mean_se_bc 
        %%Auxuliar Var names used in option B%% 'Empirical OLS Sd', 'Mean OLS Se','Mean Bootstrap Se','Bias-Corrected Rho Se',

        fprintf('======================== Summary of Results for a Sample size of %d and ρ = %.1f ========================\n',T, rho);
        disp(summ);
    
        if rho == rho_vec(1)                                                                            %Different subplots for different values of rho
            % Distribution & Bias (rho_hat vs rho_BC)
            nexttile(tl1, p);
            histogram(rho_hat_all,'Normalization','pdf','FaceAlpha',0.35,'EdgeColor','none'); hold on   %Overlaying historgrams
            histogram(rho_bc_all, 'Normalization','pdf','FaceAlpha',0.35,'EdgeColor','none');
            xline(rho,'-','LineWidth',1.5);                                                             %Adding vertival lines
            xline(mean_rho_hat,':','LineWidth',1.5);
            xline(mean_rho_bc,'-.','LineWidth',1.5);
            title(sprintf('T = %d',T))                                                                  %Adding title and labels
            xlabel('\rho estimate'); ylabel('Density');
            legend({'$\hat{\rho}$ (OLS)', '$\hat{\rho}^{BC}$', 'true $\rho$','mean $\hat{\rho}$', ...
                    'mean $\hat{\rho}^{BC}$'},'Interpreter','latex','Location','best');
            if p <= 2                                                                                   %Controlling x-axis limits
                xlim([-1.5,2])
            elseif (2 < p) & (p <= 4)
                xlim([-0.7,1.1])
            else
                xlim([-0.4,0.8])
            end
            grid on
        
            % Mean OLS SE vs Empirical SD vs Mean Bootstrap SE
            nexttile(tl2, p);
            bar_vals = [mean_se_ols, se_rho_hat, mean_se_boot_all];                                      %%Extra var option B: mean_se_bc%%
            ba = bar(categorical({'Standard Errors/Deviations'}), bar_vals);                             %Barplot of 2 se and sd
            ylabel('Standard Errors and Standard Deviation')                                             %Adding title and labels
            title(sprintf('T = %d',T))
            legend({'Mean OLS SE', 'Empirical SD($\hat{\rho}$)', 'Mean bootstrap SE'},'Interpreter', ... %%Extra var name option B:'Mean Bias Corrected $\hat{\rho}^{BC}$ SE'%%
                    'latex','Location','best');
            if p == 1                                                                                    %Controlling y-axis limits
                ylim([0,0.6])
            elseif (1 < p) & (p <= 4)
                ylim([0,0.35])
            else
                ylim([0,0.2])
            end
            grid on
        else
            % Distribution & Bias (rho_hat vs rho_BC)
            nexttile(tl3, p);
            histogram(rho_hat_all,'Normalization','pdf','FaceAlpha',0.35,'EdgeColor','none', ...
                                  'EdgeColor','none','FaceColor',[0.0 0.2 0.4]); hold on                %Overlaying historgrams
            histogram(rho_bc_all, 'Normalization','pdf','FaceAlpha',0.35,'EdgeColor','none', ...
                                  'EdgeColor','none','FaceColor',[0.6000 0.2000 0.2000]);
            xline(rho,'-','LineWidth',1.5);                                                             %Adding vertival lines
            xline(mean_rho_hat,':','LineWidth',1.5);
            xline(mean_rho_bc,'-.','LineWidth',1.5);
            title(sprintf('T = %d',T))                                                                  %Adding title and labels
            xlabel('\rho estimate'); ylabel('Density');
            legend({'$\hat{\rho}$ (OLS)', '$\hat{\rho}^{BC}$', 'true $\rho$','mean $\hat{\\rho}$', ...
                    'mean $\hat{\rho}^{BC}$'},'Interpreter','latex','Location','best');
            if p <= 2                                                                                   %Controlling x-axis limits
                xlim([-1.5,2])
            elseif (2 < p) & (p <= 4)
                xlim([0,1.2])
            else
                xlim([0.3,1])
            end
            grid on
        
            % Mean OLS SE vs Empirical SD vs Mean Bootstrap SE
            nexttile(tl4, p);
            bar_vals = [mean_se_ols, se_rho_hat, mean_se_boot_all];                                      %%Extra var option B: mean_se_bc%%
            bb = bar(categorical({'Standard Errors/Deviations'}), bar_vals);                             %Barplot of 2 se and sd
            ylabel('Standard Errors and Standard Deviation')                                             %Adding title and labels
            bb(1).FaceColor = [0.0 0.2 0.4];                                                             %Color Bar blue
            bb(2).FaceColor = [0.6000 0.2000 0.2000];                                                    %Color Bar red
            bb(3).FaceColor = [0.8500 0.6000 0.2000];                                                    %Color Bar ochre
            %bb(4).FaceColor = [0.55 0.35 0.38];                                                         %%Used in option B%% 
                                                                                                         %Color Dark Purple
            title(sprintf('T = %d',T))
            legend({'Mean OLS SE', 'Empirical SD($\hat{\rho}$)', 'Mean bootstrap SE'},'Interpreter', ... %%Extra var name option B:'Mean Bias Corrected $\hat{\rho}^{BC}$ SE'%%
                    'latex','Location','best');
            if p == 1                                                                                    %Controlling y-axis limits
                ylim([0,0.5])
            elseif (1 < p) & (p <= 4)
                ylim([0,0.3])
            else
                ylim([0,0.15])
            end
            grid on
        end
    end
end


%% OPTION A%%
%%Explanation: We consider that the bias correction is a deterministic
%%shift of the center and thys we can still use rho's varaince (ols). It is
%%simpler but comes at the cost of ignoring the extra uncertainty from
%%estimating the bias.
%%create a bootstrap analog for our BC estimation, we compute its se.

% Why is our coverage lower for BC than for OLS?
% OLS se is lower than the actual BC se, which means that our CI are
% naively narrow wich is why the coverage is also lower (more points are
% outside the bounds). 

%% OPTION B%%
%%Explanation: We manually try to use the frmula found in Efron, 1979  to
%%create a bootstrap analog for our BC estimation, we compute its se.
%%Closer to reality, se becomes very big due to the bias inclusion which is 
%%essentially a random variable,thus CI becames huge and thus more concervative
%%(reaching a coverage of 100% sometimes); on top of that, it adds complexity 
%%to the code and interpretation.

% Why are our BC se significantly higher than for OLS (why are our BC CI wider than those of OLS)?
% rho^{BC} is computed by substracting the bias=ρ^​∗​−ρ^ to rho_hat, given that
% the bias is a random variable (it depends on bootstrap draws), adds extra variability. 
% We face a clear trade-off increased variability in exchange of a less biase estimator, 
% that is why CI for the BC estimator are significantly widder.