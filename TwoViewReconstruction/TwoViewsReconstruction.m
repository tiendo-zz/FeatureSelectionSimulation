function [c2_q_c1, c2_t_c1, c1_f_hat,inlier_id] = TwoViewsReconstruction(data, max_iter)

c2_homo_full = data.C2_homo; 
c1_homo_full = data.C1_homo;
c2_q_c1 = data.C2_q_C1_init;
c2_t_c1 = data.C2_t_C1_init;


%% Function two views reconstruction:
%  Take as input the initial state estimate = [q, t, c1_p_f ...]
%  Refine the solution (including the quaternion and unit vector of translation)
%  so as to minimize the reprojection error between 2 views,
if size(c2_homo_full,2) ~= size(c1_homo_full,2)
    disp('Program terminated due to matching length is incorrect!')
    return;
end

%% GetInliers & 
N = size(c2_homo_full,2);
c1_f_hat = [];
c2_homo = [];
c1_homo = [];
error_tolerance = 0.1;
c2_R_c1 = quat2rot(c2_q_c1);
inlier_id = [];
for i = 1:N    
    c2_b_f = c2_homo_full(:,i); c2_b_f = c2_b_f / norm(c2_b_f);
    c1_b_f = c1_homo_full(:,i); c1_b_f = c1_b_f / norm(c1_b_f);
    depths = EstimateDepth(c2_b_f, c1_b_f, c2_q_c1, c2_t_c1);
    if(depths(1) > 0 && depths(2) > 0)
        c1_f_i = depths(1)*c1_b_f; 
        c2_f_i = depths(2)*c2_b_f; 
        
        c2_f_i_hat = c2_R_c1 * c1_f_i + c2_t_c1; 
        c2_f_i_hat = c2_f_i_hat / c2_f_i_hat(3);
        
        c1_f_i_hat = c2_R_c1' * (c2_f_i - c2_t_c1);
        c1_f_i_hat = c1_f_i_hat / c1_f_i_hat(3);
        err = norm(c1_f_i_hat - c1_homo_full(:,i)) + norm(c2_f_i_hat - c2_homo_full(:,i));
        if err < error_tolerance
            c1_f_hat = [c1_f_hat; depths(1)*c1_b_f];
            c2_homo = [c2_homo c2_homo_full(:,i)];
            c1_homo = [c1_homo c1_homo_full(:,i)];
            inlier_id = [inlier_id i];
        end
    end
end


N = size(c2_homo,2);

if N < 15
    disp('Not enough inliers!');
    return;
end

c2_q_c1 = data.C2_q_C1_init;
c2_t_c1 = data.C2_t_C1_init;
x_hat = [c2_q_c1; c2_t_c1];

%% Stacking features on top of each others
% for i = 1:N    
%     c2_b_f = c2_homo(:,i); c2_b_f = c2_b_f / norm(c2_b_f);
%     c1_b_f = c1_homo(:,i); c1_b_f = c1_b_f / norm(c1_b_f);
%     depths = EstimateDepth(c2_b_f, c1_b_f, c2_q_c1, c2_t_c1);
%     c1_f_hat = [c1_f_hat; depths(1)*c1_b_f];
% end


%% Reprojection error init
% [reproj_err, c2_f_hat] = Reproj_err_batch(c2_q_c1, c2_t_c1, c1_f_hat, c2_homo, c1_homo);
% prev_reproj_err = reproj_err;


%% Iterative least squares refinement
% c1_f_hat = zeros(3*N,1);
c2_f_hat = zeros(3*N,1);


H_f = zeros(4*N,3);   % Jacobian corresponding to each feature in state vector
H_r = zeros(2*N,5);   % Jacobian corresponding to the rel pose applying for each feature
H_r_small = zeros(N,5);   % Jacobian corresponding to the rel pose marginalizing features
J1_f = zeros(2*N,3);  % Image Jacobian for feature in view 1
J2_f = zeros(2*N,3);  % Image Jacobian for feature in view 2
c1_uv_tilda = zeros(2*N,1);
c2_uv_tilda = zeros(2*N,1);
r_residual = zeros(N,1);
c1_f_tilde = zeros(3*N,1);

err = 1.0;
prev_err = 0;iter = 1;
while iter < max_iter && abs(err - prev_err) > 1e-5
    [t_p, t_pp] = Orthonormal_3D_Set(c2_t_c1);
    c2_R_c1 = quat2rot(c2_q_c1);
    prev_err = err;
    
    % Loop through all features and stacking features, building Jacobians
    for i = 1:N        
        % Step 1: Stacking features
%         c2_b_f = c2_homo(:,i); c2_b_f = c2_b_f / norm(c2_b_f);
%         c1_b_f = c1_homo(:,i); c1_b_f = c1_b_f / norm(c1_b_f);
%         depths = EstimateDepth(c2_b_f, c1_b_f, c2_q_c1, c2_t_c1);
%         c1_f_hat(3*(i-1)+1:3*i) = depths(1)*c1_b_f;

        % Step 2: Calculate estimate for feature position i in camera frame 2    
        c2_f_hat(3*(i-1)+1:3*i) = c2_R_c1 * c1_f_hat(3*(i-1)+1:3*i) + c2_t_c1;

        % Step 3: Calculate measurement residual: z_tilda = meas - estimate
        c1_uv_tilda(2*(i-1)+1:2*i) = c1_homo(1:2,i) - c1_f_hat(3*(i-1)+1:3*i-1)/c1_f_hat(3*i);
        c2_uv_tilda(2*(i-1)+1:2*i) = c2_homo(1:2,i) - c2_f_hat(3*(i-1)+1:3*i-1)/c2_f_hat(3*i);        
        
        % Step 3: Build Jacobian
        J1_f(2*(i-1)+1:2*i, :) = Img_Jacobian(c1_f_hat(3*(i-1)+1:3*i));
        J2_f(2*(i-1)+1:2*i, :) = Img_Jacobian(c2_f_hat(3*(i-1)+1:3*i));
        H_r(2*(i-1)+1:2*i, :) = J2_f(2*(i-1)+1:2*i, :) * [skewsymm(c2_R_c1*c1_f_hat(3*(i-1)+1:3*i)) -t_p t_pp];
        H_f(4*(i-1)+1:4*i, :) = [J1_f(2*(i-1)+1:2*i, :); J2_f(2*(i-1)+1:2*i, :)*c2_R_c1];

        % Step 4: Left null space for the pose part first
        lns_H_f_i = zeros(4,1);
        lns_H_f_i(3) = -H_f(4*i, :)*c1_homo(:,i)/(H_f(4*i-1, :)*c1_homo(:,i));
        lns_H_f_i(4) = 1;
        lns_H_f_i(1:2) = -H_f(4*i-1:4*i, 1:2)'*lns_H_f_i(3:4);
        lns_H_f_i(3:4) = lns_H_f_i(3:4)*1/c1_f_hat(3*i);
        lns_H_f_i = lns_H_f_i/norm(lns_H_f_i);
            % Check correctness of lns:
            % lns_H_f_i'*H_f(4*(i-1)+1:4*i, :)
        
        % Step 5: Compute pose residual, pose Jacobian:
        r_residual(i) = lns_H_f_i'*[c1_uv_tilda(2*(i-1)+1:2*i);c2_uv_tilda(2*(i-1)+1:2*i)];
        H_r_small(i,:) = lns_H_f_i(3:4)'*H_r(2*(i-1)+1:2*i, :);
    end
    
    % fprintf('iter: %d, err: %f\n', iter, norm([c1_uv_tilda;c2_uv_tilda]));
    %% Update for the pose part
    [Q_r, R_r] = qr(H_r_small, 0);
    D_r = diag(R_r);
    if min(abs(D_r)) < 0.005
        disp('Cannot update new pose');
        return;
    end
    x_tilde = R_r \ (Q_r'*r_residual);
    delta_theta = [0.5*x_tilde(1:3);1]; 
    
    % Quaternion update    
    % Either normalize delta theta
    delta_theta = delta_theta / norm(delta_theta);
    c2_q_c1 = quat_mul(delta_theta, c2_q_c1);
    % Or normalize the quaternion afterward
    % x_hat(1:4) = quat_mul(delta_theta, x_hat(1:4));
    % x_hat(1:4) = x_hat(1:4)/norm(x_hat(1:4));
    
    % Position vector update
    c2_t_c1 = c2_t_c1 - x_tilde(4)*t_p + x_tilde(5)*t_pp;
    c2_t_c1 = c2_t_c1 / norm(c2_t_c1); 
    
    
    %% Solve for each features
    % Take out the pose residual part
    c2_uv_tilda = c2_uv_tilda - H_r*x_tilde;
    
    % Solve for each feature in camera 1 and 2:
    for i = 1:N
        [c1_Q_f_i, c1_R_f_i] = qr(H_f(4*(i-1)+1:4*i, :), 0);      
        test_cond = [c1_R_f_i(1,1);c1_R_f_i(2,2);c1_R_f_i(3,3)];
        if(min(abs(test_cond)) > 1e-6)
            %fprintf('update feature %d\n', i);
            c1_f_tilde(3*(i-1)+1:3*i) = c1_R_f_i \ (c1_Q_f_i'*[c1_uv_tilda(2*(i-1)+1:2*i);c2_uv_tilda(2*(i-1)+1:2*i)]);
        else
            %% add flag
            return;
%             Q_z_tilde = c1_Q_f_i'*[c1_uv_tilda(2*(i-1)+1:2*i);c2_uv_tilda(2*(i-1)+1:2*i)];
%             if (min(abs(Q_z_tilde)) < 0.0001)
%                A = [eye(3) c1_R_f_i(1:2,:)';c1_R_f_i(1:2,:) zeros(3)];
%                b = [zeros(3,1);Q_z_tilde(1:2)];
%                c1_f_tilde_lamda = A\b;
%                c1_f_tilde(3*(i-1)+1:3*i) = c1_f_tilde_lamda(1:3);
%             else
%                 return;
%             end
        end
    end
    
    % Update for the features part
    c1_f_hat = c1_f_hat + c1_f_tilde;
    
    % Update cost function value
    err = norm([c1_uv_tilda;c2_uv_tilda]);
    fprintf('iter: %d, error: %f\n', iter, err);
    iter = iter + 1;
    
%     if (prev_err < err)
%         break;
%     end
end

