% Simulation of Compute and Forward with arbitrary coarse and fine
% lattices.
% Author: Yihua Tan
% Email: adherentx.tyh@gmail.com
% The Chinese University of Hong Kong

clear all;
clc;

L = 2; % L transmitters
M = 2; % M relays
p = 7; % The prime number

% find the optimal alpha
alpha_find = @(h_m, P_mat, a_m) h_m*P_mat*P_mat'*a_m'/(1+h_m*P_mat*P_mat'*h_m');

%% Cell 1: Simulation in Only the First Hop, Equal Power Constraint
P_eq_dB_Min =40;
P_eq_dB_Max = 60;
P_eq_dB = P_eq_dB_Min:5:P_eq_dB_Max; % Power constraint in dB
P_eq = 10.^(P_eq_dB/10); % Power constraint
Pl_con = P_eq;

iter_H = 1;
sum_rate_P = zeros(length(P_eq), 1);
for i_P = 1:length(P_eq)
    for i_H = 1:iter_H
        % H_a = randn(M, L); % The channel matrix of the first hop
        H_a = [0.5 1; 2 -1];
        sum_rate_i_H = 0;
        sum_rate_i_H_var = 0;
        A_best = zeros(M, L);
        P_vec_best_for_A_best = zeros(L, 1);
        max_a = 4;
        alpha_opt_for_A_best = zeros(M, 1);
        for i_A = 1:(2*max_a)^(M*L) % try to find a good coefficient matrix A
            k = i_A;
            A_temp = zeros(M, L);
            for i_M_L = 1:M*L
                A_temp(i_M_L) = floor(k/(2*max_a)^(M*L-i_M_L));
                k = k - A_temp(i_M_L)*(2*max_a)^(M*L-i_M_L);
            end
            A = A_temp - max_a*ones(M, L);
            
            % Fixed power
            P_mat = sqrt(Pl_con(i_P))*eye(L);
            P_vec = Pl_con(i_P)*ones(L, 1);
            
            alpha_opt = zeros(M, 1);
            if sum(sum(A==[0 0; 2 -1])) == 4
            end
            for i_alpha = 1:M
                alpha_opt(i_alpha) = alpha_find(H_a(i_alpha, :), P_mat, A(i_alpha, :));
            end % for i_alpha
            
            r = zeros(L, 1);

            for i_l = 1:L
                if isempty(find(A(:, i_l) ~= 0, 1))
                    r(i_l) = 0; % all coefficients are 0.
                else
                    phi_max = 0;
                    for i_m = 1:M
                        if A(i_m, i_l) ~= 0
                            sum_mis = 0;
                            for i_mis=1:L
                                sum_mis = sum_mis+(alpha_opt(i_m)*H_a(i_m, i_mis)-A(i_m, i_mis))^2*P_vec(i_mis);
                            end % for i_mis
                            phi = alpha_opt(i_m)^2+sum_mis;
                            phi_max = max(phi, phi_max);
                        end % if A(i_m, i_l)
                    end % for i_m
                    % display(A);
                    r(i_l) = 0.5*log(max(1, P_vec(i_l)/phi_max));
                end % if isempty(find(A(:, i_l)>0)
            end % for i_l

            sum_rate_i_A = sum(r);
            if sum_rate_i_H < sum_rate_i_A
                sum_rate_i_H = sum_rate_i_A;
                A_best = A;
                alpha_opt_for_A_best = alpha_opt;
            end % if sum_rate_i_H < sum_rate_i_A
            
            % Variable power
            sum_rate_i_A_var = 0;
            P_vec_best = zeros(size(P_vec));
            division_P = 5; % divided into ? level.
            delta_P = Pl_con(i_P)/(division_P-1);
            for i_P_prod = 0:division_P^L-1
                P_prod = i_P_prod;
                for i_dim_P = 1:L
                    P_prod_t_dim = division_P^(L-i_dim_P);
                    P_temp = floor(P_prod/P_prod_t_dim)*delta_P;
                    P_mat(L-i_dim_P+1, L-i_dim_P+1) = sqrt(P_temp);
                    P_vec(L-i_dim_P+1) = P_temp;
                    P_prod = P_prod-floor(P_prod/P_prod_t_dim)*P_prod_t_dim;
                end % for i_dim_P
                alpha_opt = zeros(M, 1);
                
                for i_alpha = 1:M
                    alpha_opt(i_alpha) = alpha_find(H_a(i_alpha, :), P_mat, A(i_alpha, :));
                end % for i_alpha
                
                r = zeros(L, 1);
                for i_l = 1:L
                    if isempty(find(A(:, i_l) ~= 0, 1))
                        r(i_l) = 0; % all coefficients are 0.
                    else
                        phi_max = 0;
                        for i_m = 1:M
                            if A(i_m, i_l) ~= 0
                                sum_mis = 0;
                                for i_mis=1:L
                                    sum_mis = sum_mis+(alpha_opt(i_m)*H_a(i_m, i_mis)-A(i_m, i_mis))^2*P_vec(i_mis);
                                end % for i_mis
                                phi = alpha_opt(i_m)^2+sum_mis;
                                phi_max = max(phi, phi_max);
                            end % if A(i_m, i_l)
                        end % for i_m
                        % display(A);
                        r(i_l) = 0.5*log(max(1, P_vec(i_l)/phi_max));
                    end % if isempty(find(A(:, i_l)>0)
                end % for i_l
                
                if sum_rate_i_A_var < sum(r)
                    sum_rate_i_A_var = sum(r);
                    P_vec_best = P_vec;
                end
            end % for i_P_prod
%             if sum_rate_i_H_var < sum_rate_i_A_var
%                 sum_rate_i_H_var = sum_rate_i_A_var;
%                 A_best = A;
%                 P_vec_best_for_A_best = P_vec_best;
%             end % if sum_rate_i_H_var < sum_rate_i_A_var

            if false
            % Compare variable and fixed power methods when A is fixed.
            if sum_rate_i_A_var_P_mat > sum_rate_i_A
                display(['Variable power is better. ', ...
                    'sum_rate_i_A_var=' num2str(sum_rate_i_A_var), ' ', ...
                    'sum_rate_i_A=' num2str(sum_rate_i_A), ' ']);
                display(P_vec_best);
%                 display(['Power Constraint: ', num2str(Pl_con(i_P))]);
%                 display(A);
%                 display(H_a);
            elseif sum_rate_i_A_var_P_mat < sum_rate_i_A
                display(['Fixe power is better. ', ...
                    'sum_rate_i_A_var=' num2str(sum_rate_i_A_var), ' ', ...
                    'sum_rate_i_A=' num2str(sum_rate_i_A), ' ']);
                display(A);
                display(H_a);
                error('This is an unexpected result.');
            end
            end % if commented
        end % for i_A
        
        display(H_a);
        display(alpha_opt_for_A_best);
        display(A_best);
        
        % Compare variable and fixed power methods when H is fixed.
        if false
        if sum_rate_i_H_var > sum_rate_i_H
            display(['Variable power is better. ', ...
                'sum_rate_i_H_var=' num2str(sum_rate_i_H_var), ' ', ...
                'sum_rate_i_H=' num2str(sum_rate_i_H), ' ']);
            display(H_a);
            display(A_best);
            display(P_vec_best_for_A_best);
        elseif sum_rate_i_H_var < sum_rate_i_H
            display(['Fixed power is better. ', ...
                'sum_rate_i_H_var=' num2str(sum_rate_i_H_var), ' ', ...
                'sum_rate_i_H=' num2str(sum_rate_i_H), ' ']);
            display(H_a);
            display(A_best);
            display(P_vec_best_for_A_best);
            error('This is an unexpected result.');
        end % if sum_rate_i_H_var > sum_rate_i_H
        end % if commented
    end % for i_H
    sum_rate_P(i_P) = 0;
end % for i_P
%plot(P_eq_dB, sum_rate_P);






