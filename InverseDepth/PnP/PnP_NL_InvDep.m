function [Ci_R_Cr, Ci_p_Cr, normHist] = PnP_NL_InvDep(Ci_R_Cr, Ci_p_Cr, AbsolutePoses, Cr_p_f, Ci_z_f, CameraParams)

maxIter = 100;

Ci_z_f_norm = CameraParams.Kinv*[Ci_z_f; ones(1,size(Ci_z_f,2))];
Ci_z_f_norm = Ci_z_f_norm(1:2,:) ./ repmat(Ci_z_f_norm(3,:), [2 1]);

% Gauss-Newton iteration
Cr_p_Ci_iter = -Ci_R_Cr'*Ci_p_Cr;
Ci_q_Cr_iter = rot2quat(Ci_R_Cr);
iter = 1;
norm_res_log = zeros(1,maxIter);
res = 1;
while norm(res) > 1e-5 && iter < maxIter
    Jpose = zeros(2*size(Cr_p_f,2),6);
    for featk = 1:size(Cr_p_f,2)
        Cr_p_f_single = Cr_p_f(1:3 ,featk);              
   
        % Pose Transformation
        Cr_T_W = AbsolutePoses(:,:,Cr_p_f(4 ,featk));
        W_R_Cr = Cr_T_W(:,1:3)';
        W_p_Cr = -Cr_T_W(:,1:3)'*Cr_T_W(:,4);
        
        [Jq,res_tmp] = dprojection_dq(Cr_p_f_single, quat2rot(Ci_q_Cr_iter), W_R_Cr,W_p_Cr, Cr_p_Ci_iter, Ci_z_f_norm(1:2,featk));
        Jp = dprojection_dp(Cr_p_f_single, quat2rot(Ci_q_Cr_iter), W_R_Cr,W_p_Cr, Cr_p_Ci_iter, Ci_z_f_norm(1:2,featk));
        Jpose(2*featk-1:2*featk,:) = [Jp Jq]; 
        res(2*featk-1:2*featk,:) = res_tmp;
    end
    
        
    
    % Solver:
    dx = (Jpose'*Jpose + 1*diag([0.1 0.1 0.1 1 1 1])) \ (Jpose'*res);
    
    step_size = min(1, 0.9/norm(dx(4:6)));
    dx = dx*step_size;
    dp = dx(1:3); dtt = dx(4:6);
    
    % update:
    Cr_p_Ci_iter = Cr_p_Ci_iter + dp;
    dq = [1/2*dtt; sqrt(1 - 1/4*norm(dtt)^2)];
    Ci_q_Cr_iter = quat_mul(dq, Ci_q_Cr_iter);
    
    
    % Logger
    norm_res_log(iter) = norm(res);
    iter = iter + 1;
    
%    fprintf('Iter :%d, norm: %f\n', iter, norm(res));

end

Ci_R_Cr = quat2rot(Ci_q_Cr_iter);
Ci_p_Cr = -Ci_R_Cr*Cr_p_Ci_iter;


normHist = norm_res_log(1:iter-1);