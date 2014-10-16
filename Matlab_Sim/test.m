            for i_l = 1:L
                a_xl = A(:, i_l);
                idx_pos_a = find(a_xl ~= 0);
                if isempty(idx_pos_a)
                    r(i_l) = 0;
                else
                    sum_mis = zeros(M, 1);
                    for i_mis=1:L
                        sum_mis = sum_mis+(alpha_opt.*H_a(:, i_mis)-A(:, i_mis)).^2*P_vec(i_mis);
                    end % for i_mis
                    phi = alpha_opt.^2+sum_mis;
                    phi_max = max(phi(idx_pos_a));
                    r(i_l) = 0.5*log(max(1, P_vec(i_l)/phi_max));
                end % if isempty(idx_pos_a)
            end % for i_l