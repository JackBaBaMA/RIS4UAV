clc
% variables declaration
alpha = -43;
sigma_2 = -70;
lambda_wave = 0.087;
l_m = lambda_wave / 10;
l_n = lambda_wave / 2;
d_move = 10;    % moving distance in one slot

gamma_1 = 1;
gamma_2 = 1;
gamma1_star = 1;
gamma2_star = 1;

v_1a = [-1/2, sqrt(3)/2,0];   % base direction of RIS 1
v_1b = [-1/2, sqrt(3)/2,0];   % base direction of RIS 2
n_r = 2; % antena number of transitter
n_t = 2; % antena number of receiver

% Assuming u_i ,u_r and u_t are vectors (1D NumPy arrays)
u_1 = [0, 0, 0];
u_2 = [0, 50, 0]; % location of UAV
u_t = [1, 0, 0];  % Replace u_t_x, u_t_y, u_t_z  with actual values
u_r = [1, 50, 0];
v_t = [-1, 0, 0];
v_r = [-1, 0, 0];

% Calculate the Euclidean distance
d_T1 = norm(u_1 - u_t);
d_T2 = norm(u_2 - u_t);
d_R1 = norm(u_1 - u_r);
d_R2 = norm(u_2 - u_r);
d_S = norm(u_1 - u_2);
achievable_rate_p = [];

p = 5; % power

%for element_num = 2:2:60
element_num = 2;
    for n = 1:10000
        if n == 1
            u_2_old = u_2;
        else
            u_2_old = u_2_optimal;
        end

        u_0 = rand();    % initialize   u_0, y_0 and w_0
        v_0 = rand();
        w_0 = rand();




        cvx_begin
            variable u_2(3);     % optimize the location of UAV, a vector with 3 elements
            variable u
            variable v
            variable w
    
            % step 1: get Q *
            M_1 = element_num;
            M_2 = element_num;
            [r_1L_list, r_2L_list, t_1R_list, t_2R_list] = get_r1l(u_1, u_r, v_r, u_2, u_t, v_t, lambda_wave, n_r, n_t);
            h_n = get_h(alpha, M_1, M_2, d_R1, d_T1, d_R2, d_T2, d_S, lambda_wave, gamma_1, gamma_2, r_1L_list, r_2L_list, t_1R_list, t_2R_list);
            [U, S, V] = svd(h_n);
            delta = diag(S);
            Q_star = compute_optimal_power_allocation(sigma_2, delta, p, V);
    
            % step 2: get gamma 1
            [U, Sigma, V] = svd(Q_star);
            r = rank(Q_star);
            U_q = U(:, 1:r);
            Sigma_q = diag(Sigma(1:r));
    
            beta_a = exp(-1j * 2 * pi * (d_R1 + d_T2) / lambda_wave);
            beta_b = exp(-1j * 2 * pi * (d_R2 + d_T2) / lambda_wave);
            beta_c = exp(-1j * 2 * pi * (d_R2 + d_S + d_T1) / lambda_wave);
    
            A_n = U_q * sqrt(Sigma_q);
            B_n = (alpha * M_2 * beta_b) / (d_R2 * d_T2) * (r_2L_list.' * t_2R_list);
            C_n = ((alpha^(3/2)) * M_1 * M_2 * beta_c) / (d_R2 * d_S * d_T1) * (r_2L_list.' * t_1R_list);
    
            x1_n = ones(1, n_r) + (1 / sigma_2) * (A_n * A_n' + B_n * B_n' + C_n * C_n' + gamma_2 * C_n * A_n' + gamma_2 * A_n * C_n');
            y1_n = (1 / sigma_2) * (gamma_2 * A_n * B_n' + gamma2_star * C_n * B_n');
            x1_n_inv = inv(x1_n);
            v1 = x1_n_inv * y1_n;
            gamma1_star = getgamma_star(y1_n, x1_n_inv, v1);
    
            % step 3: get gamma2*
            x2_n = ones(1, n_r) + (1 / sigma_2) * (A_n * A_n' + B_n * B_n' + C_n * C_n' + gamma_1 * C_n * B_n' + gamma1_star * B_n * C_n');
            y2_n = (1 / sigma_2) * (gamma1_star * B_n * A_n' + C_n * A_n');
            x2_n_inv = inv(x2_n);
            v2 = x2_n_inv * y2_n;
            gamma2_star = getgamma_star(y2_n, x2_n_inv, v2);
    
            % step 4: get L
            x3_n = alpha * M_1 * beta_a * gamma1_star * (r_1L_list.' * t_1R_list);
            y3_n = alpha * M_2 * beta_b * gamma2_star * (r_2L_list.' * t_2R_list);
            z_n = (((alpha^(3/2)) * M_1 * M_2 * beta_c * gamma1_star * gamma2_star) / d_T1) * (r_2L_list.' * t_1R_list);
   
    
            A_0_inside = 1 + (1 / sigma_2) * (x3_n^2 + y3_n^2 / (u * v) + z_n^2 / (u * w) + ...
                2 * x3_n * y3_n / (sqrt(u) * sqrt(v)) + ...
                2 * x3_n * z_n / (sqrt(u) * sqrt(w)) + ...
                2 * y3_n * z_n / (u * sqrt(v) * sqrt(w)));
    
            A_0 = log2(A_0_inside);
    
            B_0 = -(1 / sigma_2) * (y3_n^2 / (u^2 * v) + z_n^2 / (u^2 * w) + ...
                x3_n * z_n / (sqrt(u)^3 * sqrt(w)) + ...
                4 * y3_n * z_n / (u^2 * sqrt(v) * sqrt(w)));
    
            C_0 = -(1 / sigma_2) * (y3_n^2 / (u * v^2) + x3_n * y3_n / (sqrt(u) * sqrt(v)^3) + ...
                y3_n * z_n / (u * sqrt(v)^3 * sqrt(w)));
    
            D_0 = -(1 / sigma_2) * (z_n^2 / (u * w^2) + x3_n * z_n / (sqrt(u) * sqrt(w)^3) + ...
                y3_n * z_n / (u * sqrt(v) * sqrt(w)^3));
    
            obj = B_0 / (A_0 * log(2)) * u + C_0 / (A_0 * log(2)) * v + D_0 / (A_0 * log(2)) * w;
            maximize(obj)
            subject to
                quad_form(u) <= d_R2^2 + u_0^2 - 2 * u_0 * u;
                quad_form(v) <= d_T2^2 + v_0^2 - 2 * v_0 * v;
                quad_form(w) <= d_S^2 + w_0^2 - 2 * w_0 * w;
                norm(u_2 - u_2_old) <= d_move;




        % Use u_2_optimal as needed
        disp(['Optimal value of u_2: ', num2str(u_2_optimal)]);


        % Calculate the trace of the matrix expression
        trace_expr = trace(h_n * Q_star * h_n');

        % Define the objective function
        achievable_rate = log2(det(ones(n_r) + (1 / sigma_2) * h_n * Q_star * h_n')).real;
        achievable_rate_p = [achievable_rate_p, abs(achievable_rate)];

    end
%end






