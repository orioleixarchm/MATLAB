% Filename: sticky_prices_wages_loop_matlab.m

% List of iterative variables
Theta_p   = [0.75; 0.25; 0.75];
Theta_w   = [0.75; 0.75; 0.25];
N       = length(Theta_p);
save pars1 Theta_p Theta_w N

% Creating IRFs 
IRFs = zeros(20,N,9);
save IRF_for_looping IRFs

% Looping and solving the model for each parameter value
for j=1:N
    load pars1
    curr_theta_p = Theta_p(j);
    curr_theta_w = Theta_w(j);
    save pars2 curr_theta_p curr_theta_w j
    dynare sticky_prices_wages_loop;
    load IRF_for_looping
    load pars2
    IRFs(:,j,1) = oo_.irfs.y_tilde_epsilon_a;
    IRFs(:,j,2) = oo_.irfs.y_tilde_epsilon_z;
    IRFs(:,j,3) = oo_.irfs.y_tilde_epsilon_nu;
    IRFs(:,j,4) = oo_.irfs.pi_p_epsilon_a;
    IRFs(:,j,5) = oo_.irfs.pi_p_epsilon_z;
    IRFs(:,j,6) = oo_.irfs.pi_p_epsilon_nu;
    IRFs(:,j,7) = oo_.irfs.pi_w_epsilon_a;
    IRFs(:,j,8) = oo_.irfs.pi_w_epsilon_z;
    IRFs(:,j,9) = oo_.irfs.pi_w_epsilon_nu;
    save IRF_for_looping IRFs
end

% visualization
for j = 1:N
    leg{j} = sprintf('theta_p %.2f, theta_w %.2f', Theta_p(j), Theta_w(j)); 
end

figure
subplot(3,3,1)
plot(IRFs(:,:,1))
title('y\_tilde\_epsilon\_a')
legend(leg,'Location','SouthEast')
subplot(3,3,2)
plot(IRFs(:,:,2))
title('y\_tilde\_epsilon\_z')
legend(leg,'Location','SouthEast')
subplot(3,3,3)
plot(IRFs(:,:,3))
title('y\_tilde\_epsilon\_nu')
legend(leg,'Location','SouthEast')
subplot(3,3,4)
plot(IRFs(:,:,4))
title('pi\_p\_epsilon\_a')
legend(leg,'Location','SouthEast')
subplot(3,3,5)
plot(IRFs(:,:,5))
title('pi\_p\_epsilon\_z')
legend(leg,'Location','SouthEast')
subplot(3,3,6)
plot(IRFs(:,:,6))
title('pi\_p\_epsilon\_nu')
legend(leg,'Location','SouthEast')
subplot(3,3,7)
plot(IRFs(:,:,7))
title('pi\_w\_epsilon\_a')
legend(leg,'Location','SouthEast')
subplot(3,3,8)
plot(IRFs(:,:,8))
title('pi\_w\_epsilon\_z')
legend(leg,'Location','SouthEast')
subplot(3,3,9)
plot(IRFs(:,:,9))
title('pi\_w\_epsilon\_nu')
legend(leg,'Location','SouthEast')