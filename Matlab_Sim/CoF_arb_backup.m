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
P_eq_dB_Max = 20;
P_eq_dB = 5:5:P_eq_dB_Max; % Power constraint in dB
P_eq = 10.^(P_eq_dB/10); % Power constraint
Pl_con = P_eq;

iter_H = 10;
iter_A = 100;
sum_rate_P = zeros(length(P_eq), 1);
for i_P = 1:length(P_eq)
    for i_H = 1:iter_H
        H_a = randn(L, M); % The channel matrix of the first hop
        sum_rate = 0;
        sum_rate_var = 0;
        for i_A = 1:iter_A % try to find a good coefficient matrix A
            % TODO: Use a better way to generate matrix A
            while true
                A = randi(p, L, M)-1; % Tentatively we randomize A.
                if isempty(find(sum(A, 2)==0, 1))
                    break;
                end
            end
            
            % TODO: Loop: optimize over P_mat
            
            % Fixed power
            P_mat = sqrt(Pl_con(i_P))*eye(L);
            P_vec = Pl_con(i_P)*ones(L, 1);
            
            alpha_opt = zeros(M, 1);
            for i_alpha = 1:M
                alpha_opt(i_alpha) = alpha_find(H_a(i_alpha, :), P_mat, A(i_alpha, :));
            end % for i_alpha
            
            r = zeros(L, 1);
            for i_r = 1:L
                a_xl = A(:, i_r);
                idx_pos_a = find(a_xl>0);
                sum_mis = zeros(M, 1);
                for i_mis=1:L
                    sum_mis = sum_mis+(alpha_opt.*H_a(:, i_mis)-A(:, i_mis)).^2*P_vec(i_mis);
                end % for i_mis
                phi = alpha_opt.^2+sum_mis;
                phi_max = max(phi(idx_pos_a));
                r(i_r) = 0.5*log(max(1, P_vec(i_r)/phi(i_r)));
            end % for i_r
            sum_rate = max(sum_rate, sum(r));
            
            % Variable power
            sum_rate_var_P_mat = 0;
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
                for i_r = 1:L
                    a_xl = A(:, i_r);
                    idx_pos_a = find(a_xl>0);
                    sum_mis = zeros(M, 1);
                    for i_mis=1:L
                        sum_mis = sum_mis+(alpha_opt.*H_a(:, i_mis)-A(:, i_mis)).^2*P_vec(i_mis);
                    end % for i_mis
                    phi = alpha_opt.^2+sum_mis;
                    phi_max = max(phi(idx_pos_a));
                    r(i_r) = 0.5*log(max(1, P_vec(i_r)/phi(i_r)));
                end % for i_r
                if sum_rate_var_P_mat < sum(r)
                    sum_rate_var_P_mat = sum(r);
                    P_vec_best = P_vec;
                end
            end % for i_P_prod
            sum_rate_var = max(sum_rate_var, sum_rate_var_P_mat);
            if sum_rate_var > sum_rate
                display(['Variable power is better. ', ...
                    'sum_rate_var=' num2str(sum_rate_var), ' ', ...
                    'sum_rate=' num2str(sum_rate), ' ']);
                display(P_vec_best);
                display(['Power Constraint: ', num2str(Pl_con(i_P))]);
                display(A);
                display(H_a);
            end
            if sum_rate_var < sum_rate
                display(['Fixe power is better. ', ...
                    'sum_rate_var=' num2str(sum_rate_var), ' ', ...
                    'sum_rate=' num2str(sum_rate), ' ']);
                display(A);
                display(H_a);
                error('This is an unexpected result.');
            end
        end % for i_A
    end % for i_H
    sum_rate_P(i_P) = sum_rate;
end % for i_P
%plot(P_eq_dB, sum_rate_P);






