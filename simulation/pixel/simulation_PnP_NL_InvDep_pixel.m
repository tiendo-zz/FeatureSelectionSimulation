function [Ci_R_W, Ci_p_W, CameraParams, normHist] = simulation_PnP_NL_InvDep_pixel(Ci_R_W, Ci_p_W, AbsolutePoses, Cr_p_f, Ci_z_f, CameraParams)

maxIter = 10;

Ci_z_f_norm = [Ci_z_f; ones(1,size(Ci_z_f,2))];
Ci_z_f_norm = Ci_z_f_norm(1:2,:) ./ repmat(Ci_z_f_norm(3,:), [2 1]);

% Gauss-Newton iteration
W_p_Ci_iter = -Ci_R_W'*Ci_p_W;
Ci_q_W_iter = rot2quat(Ci_R_W);
iter = 1;
norm_res_log = zeros(1,maxIter);
res = 10;

regulator = 1;

while norm(res) > 1 && iter < maxIter
    Jpose = zeros(2*size(Cr_p_f,2),6);
    Jfc_cc = zeros(2*size(Cr_p_f,2),4);
    for featk = 1:size(Cr_p_f,2)
        Cr_p_f_single = Cr_p_f(1:3 ,featk);              
   
        % Pose Transformation
        Cr_T_W = AbsolutePoses(:,:,Cr_p_f(4 ,featk));
        W_R_Cr = Cr_T_W(:,1:3)';
        W_p_Cr = -Cr_T_W(:,1:3)'*Cr_T_W(:,4);

        
        m_theta_phi = [cos(Cr_p_f(2,featk))*cos(Cr_p_f(1,featk));cos(Cr_p_f(2,featk))*sin(Cr_p_f(1,featk));sin(Cr_p_f(2,featk))];
        Ci_p_f = quat2rot(Ci_q_W_iter)*W_R_Cr*m_theta_phi + Cr_p_f(3,featk)*quat2rot(Ci_q_W_iter)*(W_p_Cr-W_p_Ci_iter);
        H_pixel = CameraParams.K(1:2,1:2);         
        [Jq,~] = dprojection_dq(Cr_p_f_single, quat2rot(Ci_q_W_iter), W_R_Cr,W_p_Cr, W_p_Ci_iter, Ci_z_f_norm(1:2,featk));
        Jp = dprojection_dp(Cr_p_f_single, quat2rot(Ci_q_W_iter), W_R_Cr,W_p_Cr, W_p_Ci_iter, Ci_z_f_norm(1:2,featk));
        Jpose(2*featk-1:2*featk,:) = H_pixel*[Jp Jq]; 
        
        [Jfc,Jcc,~,~] = dprojection_dfc_dcc_simulation(Cr_p_f_single, quat2rot(Ci_q_W_iter), W_R_Cr,W_p_Cr, W_p_Ci_iter, Ci_z_f_norm(1:2,featk),CameraParams,30);
        Jfc_cc(2*featk-1:2*featk,:) = [Jfc Jcc];
        % residual
        Ck_z_f_hat = CameraParams.K*(Ci_p_f / Ci_p_f(3));
        res_tmp = Ci_z_f_norm(1:2,featk) - Ck_z_f_hat(1:2);
        res(2*featk-1:2*featk,:) = res_tmp;
    end
    
    Jpose = [Jpose Jfc_cc];    
    
    % Solver:
%     dx = (Jpose'*Jpose + 1*diag([0.1 0.1 0.1 1 1 1])) \ (Jpose'*res);
    dx = (Jpose'*Jpose + regulator*diag([0.1 0.1 0.1 1 1 1 1 1 1 1])) \ (Jpose'*res);
    
    step_size = min(1, 0.9/norm(dx(4:6)));
    dx = dx*step_size;
    dp = dx(1:3); dtt = dx(4:6);
    
    % update:
    W_p_Ci_iter = W_p_Ci_iter + dp;
    dq = [1/2*dtt; sqrt(1 - 1/4*norm(dtt)^2)];
    Ci_q_W_iter = quat_mul(dq, Ci_q_W_iter);
    
    CameraParams.fc(1) = CameraParams.fc(1) + dx(1);
    CameraParams.fc(2) = CameraParams.fc(2) + dx(2);
    CameraParams.cc(1) = CameraParams.cc(1) + dx(3);
    CameraParams.cc(2) = CameraParams.cc(2) + dx(4);
    fx = CameraParams.fc(1); fy = CameraParams.fc(2);
    cx = CameraParams.cc(1); cy = CameraParams.cc(2);
    CameraParams.K = [fx 0 cx; ...
                      0 fy cy; ...
                      0 0 1];
    CameraParams.Kinv = [1/fx 0 -cx/fx; ...
                         0 1/fy -cy/fy; ...
                         0 0 1]; 
    % Logger
    norm_res_log(iter) = norm(res);
    
    
    fprintf('Iter :%d, norm: %f, regulator: %f\n', iter, norm(res),regulator);

    if iter >= 2
       if norm_res_log(iter-1) > norm_res_log(iter)
          regulator = regulator * 0.8;
       else
          regulator = regulator * 1.25;
       end
    end
    iter = iter + 1;
end


Ci_R_W = quat2rot(Ci_q_W_iter);
Ci_p_W = -Ci_R_W*W_p_Ci_iter;


normHist = norm_res_log(1:iter-1);

end