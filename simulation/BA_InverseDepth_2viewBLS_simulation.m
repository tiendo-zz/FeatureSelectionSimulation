function [FeatureBag, AbsolutePoses, normHist, gradnormHist,EraseId,error] = BA_InverseDepth_2viewBLS_simulation(FeatureBag, AbsolutePoses, PoseGraphMatrix, featureExtracted, usedPoses, CameraParams)

PoseGraphMatrix_r = PoseGraphMatrix(:, usedPoses);
NumOfPoses = length(usedPoses);
% NumOfFeatures = length(FeatureBag(FeatureBag(4,:)~=0));

% figure(3); clf; hold on;
% Draw3DWorld(AbsolutePoses, FeatureBag, [eye(3) zeros(3,1)], usedPoses);
% axis([-2 2 -2 2 -1 5]);
%% Bundle Adjustment
regulator = 1;
norm_iter = 100;
grad_norm = 100;
normHist = [];
gradnormHist = [];
maxIter = 50;
iter = 1;

    
% 1. Build Jacobian Jp Jf corresponds to each feature
featureAvailableList_all = find(FeatureBag(4,:));
poses_see_featk_all = length(find(PoseGraphMatrix_r(featureAvailableList_all, :)));
% Jpose = zeros(2*nnz(PoseGraphMatrix_r), 6*NumOfPoses);
% Jfeat = zeros(2*nnz(PoseGraphMatrix_r), 3*NumOfFeatures);
% residual = zeros(2*nnz(PoseGraphMatrix_r), 1);

featureAvailableList = find(FeatureBag(4,:));
featureAvailableList = featureAvailableList-1;


ErasePose = cell(featureAvailableList_all(end),1);
for e = 1:length(ErasePose)
    ErasePose{e} = [];
end


while iter < maxIter && norm_iter > 1e-4 && grad_norm > 1e-6
    
AbsolutePoses_tmp = AbsolutePoses;
FeatureBag_tmp = FeatureBag;
   
NumOfFeatures = length(featureAvailableList);    
poses_see_featk_all = length(find(PoseGraphMatrix_r(featureAvailableList+1, :)));


feat_cam = zeros(poses_see_featk_all,2);

Jpose = zeros(2*poses_see_featk_all, 6*(NumOfPoses-1)+5);
Jfeat = zeros(2*poses_see_featk_all, 3*NumOfFeatures);
error = zeros(poses_see_featk_all,2);
residual = zeros(2*poses_see_featk_all, 1);
meascount = 0;
last_meas_count = 0;

featcount = 0;
last_feat_count = 0;
% Hessian build up
A22_inv = zeros(3*NumOfFeatures);


EraseId = zeros(length(featureAvailableList),1);
outlier_num = 0;
for featk = featureAvailableList %0:size(PoseGraphMatrix,1)-1
   
    
   poses_see_featk = find(PoseGraphMatrix_r(featk+1, :));
   
   count_measure_single = 0;
   for l = 0:length(poses_see_featk)-1
       
      %  Pose Jacobian      
      Cr_p_f = FeatureBag(1:3 ,featk+1);      
      Ck_z_f = CameraParams.Kinv(1:2,:)*[featureExtracted{usedPoses(poses_see_featk(l+1))}(1:2,PoseGraphMatrix_r(featk+1,poses_see_featk(l+1)));1];                  
      
      
      % Pose Transformation
      Ck_T_W = AbsolutePoses(:,:,usedPoses(poses_see_featk(l+1)));
      Cr_T_W = AbsolutePoses(:,:,FeatureBag(4, featk+1));      
      Ck_T_Cr = Ck_T_W*[InversePose(Cr_T_W);zeros(1,3) 1];            
      
      Ck_R_W = Ck_T_W(:,1:3);
      W_R_Cr = Cr_T_W(:,1:3)';
      W_p_Cr = -Cr_T_W(:,1:3)'*Cr_T_W(:,4);
      W_p_Ck = -Ck_T_W(:,1:3)'*Ck_T_W(:,4);
      
      
      m_theta_phi = [cos(Cr_p_f(2))*cos(Cr_p_f(1));cos(Cr_p_f(2))*sin(Cr_p_f(1));sin(Cr_p_f(2))];
      Ck_z_f_hat = m_theta_phi(1:2) / m_theta_phi(3);
      
      ref_idx_usePoses = find(FeatureBag(4, featk+1) == usedPoses);
      
      if usedPoses(poses_see_featk(l+1))==FeatureBag(4, featk+1)
         if usedPoses(poses_see_featk(l+1)) == 2
             Jp = zeros(2,2);
         else
             Jp = zeros(2,3);             
         end
         if usedPoses(ref_idx_usePoses) == 2
             Jpr = zeros(2,2);
         else
             Jpr = zeros(2,3);
         end
         Jq = zeros(2,3);
         Jqr = zeros(2,3);
         Jro = zeros(2,1);
         residual_tmp =  Ck_z_f - Ck_z_f_hat;         
      else
          if usedPoses(poses_see_featk(l+1)) == 2
             [Jp,~,~,t_p, t_pp] = dprojection_dp_2viewBLS(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
          else
              Jp = dprojection_dp(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);              
          end
          if usedPoses(ref_idx_usePoses) == 2
              Jpr = dprojection_dpr_2viewBLS(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
          else
              Jpr = dprojection_dpr(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
          end
          [Jq,residual_tmp,~] = dprojection_dq(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
          Jqr = dprojection_dqr(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
          Jro = dprojection_dro(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
      end
      Jtheta = dprojection_dtheta(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
      Jphi = dprojection_dphi(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
      
%       
%       z_ = Ck_z_f_hat + Jq*pertube.q_tilde(:,usedPoses(poses_see_featk(l+1)))+Jp*pertube.p_tilde(:,usedPoses(poses_see_featk(l+1)))+...
%             [Jtheta Jphi Jro]*pertube.F_tilde(:,featk+1);
%       
%       error(meascount+1,1) = norm(z_ - Ck_z_f);
%       error(meascount+1,2) = norm(residual_tmp);
      
%       gamma = residual_tmp'*10^4*residual_tmp;
%       if  gamma > 16.27  % chi-square for d.o.f = 3, p = 0.001
%         continue;
%       end
      if usedPoses(poses_see_featk(l+1)) == 2
         Jpose(2*meascount+1:2*meascount+2, 6*(poses_see_featk(l+1)-1)+1:6*(poses_see_featk(l+1)-1)+5) = [Jp Jq];
      else
          if usedPoses(poses_see_featk(l+1)) ~= 1      
              Jpose(2*meascount+1:2*meascount+2, 6*(poses_see_featk(l+1)-2)+5+1:6*(poses_see_featk(l+1)-2)+5+6) = [Jp Jq];
          end
      end
      
      
      if ref_idx_usePoses ~= 1      
         if usedPoses(ref_idx_usePoses) == 2
            Jpose(2*meascount+1:2*meascount+2, 6*(ref_idx_usePoses-1)+1:6*(ref_idx_usePoses-1)+5) = [Jpr Jqr];
         else
             if usedPoses(ref_idx_usePoses) ~= 1      
                Jpose(2*meascount+1:2*meascount+2, 6*(ref_idx_usePoses-2)+5+1:6*(ref_idx_usePoses-2)+5+6) = [Jpr Jqr];
             end
         end
      end
      
      residual(2*meascount+1:2*meascount+2) = residual_tmp;
      Jfeat(2*meascount+1:2*meascount+2, 3*featcount+1:3*featcount+3) = [Jtheta Jphi Jro];              
      
      feat_cam(meascount+1,:) = [featk+1 poses_see_featk(l+1)];
      meascount = meascount + 1;
      count_measure_single = count_measure_single + 1;
   end   

   
   Jfeatk = Jfeat(last_meas_count+1:2*meascount, 3*featcount+1:3*featcount+3);
   A22_inv(3*featcount+1:3*featcount+3, 3*featcount+1:3*featcount+3) = (Jfeatk'*Jfeatk + regulator*eye(3)) \ eye(3);
      
   featcount = featcount + 1;
   last_meas_count = 2*meascount;
end

Jpose(:,1:6) = [];


% clear the all zero lines
Jfeat(:,3*featcount+1:end) = [];
Jfeat(2*meascount+1:end,:) = [];
Jpose(2*meascount+1:end,:) = [];
residual(2*meascount+1:end,:) = [];
A22_inv(3*featcount+1:end,:) = [];
A22_inv(:,3*featcount+1:end) = [];

if iter >= 2
    if norm_iter > norm(residual)
        regulator = regulator * 0.8;
    else
        regulator = regulator * 1.25;
    end
end

norm_iter = norm(residual);
% if ~isempty(normHist)
%     if norm_iter > normHist(end)
%        break;
%     end
% end
normHist = [normHist norm_iter];
% if iter > 20
%     norm_iter = min(norm_iter,1e2*abs(normHist(end-1) - normHist(end)));
% end

Jpose_sp = sparse(Jpose);
Jfeat_sp = sparse(Jfeat);
A11 = Jpose_sp'*Jpose_sp;
A12 = Jpose_sp'*Jfeat_sp;
A22_inv_sp = sparse(A22_inv);
res_ = [Jpose_sp'*residual;Jfeat_sp'*residual];

% solve:
p_tilde = (A11 +regulator*speye(size(A11))- A12*A22_inv_sp*A12') \ (res_(1:6*(NumOfPoses-2)+5) - A12*A22_inv_sp*res_(1+6*(NumOfPoses-2)+5:end));

% do not update first 3 camera poses (if more than 3 poses)
for l = 0:NumOfPoses-2
   if l == 0
       dq = [0.5*p_tilde(6*l+3:6*l+5);1];
       dq = dq / norm(dq);
       c1_t_c2 = -AbsolutePoses_tmp(:,1:3,usedPoses(l+2))'*AbsolutePoses_tmp(:,4,usedPoses(l+2));
       AbsolutePoses_tmp(:,1:3,usedPoses(l+2)) = quat2rot( quat_mul(dq, rot2quat(AbsolutePoses_tmp(:,1:3,usedPoses(l+2)))) );      
       c1_t_c2 = c1_t_c2 - p_tilde(6*l+1)*t_p + p_tilde(6*l+2)*t_pp;
       c1_t_c2 = c1_t_c2 / norm(c1_t_c2);
       AbsolutePoses_tmp(:,4,usedPoses(l+2)) = -AbsolutePoses_tmp(:,1:3,usedPoses(l+2))*c1_t_c2;
   else 
       dq = [0.5*p_tilde(6*(l-1)+5+4:6*(l-1)+5+6);1];
       dq = dq / norm(dq);
       W_p_Ck = -AbsolutePoses_tmp(:,1:3,usedPoses(l+2))'*AbsolutePoses_tmp(:,4,usedPoses(l+2)) + p_tilde(6*(l-1)+5+1:6*(l-1)+5+3);
       AbsolutePoses_tmp(:,1:3,usedPoses(l+2)) = quat2rot( quat_mul(dq, rot2quat(AbsolutePoses_tmp(:,1:3,usedPoses(l+2)))) );       
       AbsolutePoses_tmp(:,4,usedPoses(l+2)) = -AbsolutePoses_tmp(:,1:3,usedPoses(l+2))*W_p_Ck;
   end
end


grad_norm = norm(p_tilde)^2;
 
% 4. Efficient solve for each feature f
bf_p = res_(1+6*(NumOfPoses-2)+5:end) - A12'*p_tilde;
f_tilde = A22_inv_sp * bf_p;
f_tilde_old = f_tilde;
% [~,erase_b] = ismember(EraseId,featureAvailableList);
% FeatureBag(1:3, featureAvailableList(erase_b)+1) = 0; % erase outliers
% featureAvailableList(erase_b) = [];
FeatureBag_tmp(1:3, featureAvailableList+1) = FeatureBag_tmp(1:3, featureAvailableList+1) + reshape(f_tilde, [3, length(featureAvailableList)]);

%%after updating the temp result, compute the cost function one more time,
%%if it's smaller than the previous, take it and decrease regulator, else,
%%discard it and increase twice the regulator

cost_function_value = ComputeCostFunction(FeatureBag_tmp, AbsolutePoses_tmp, PoseGraphMatrix_r, featureExtracted, featureAvailableList,usedPoses,CameraParams);




grad_norm = norm(f_tilde)^2 + norm(p_tilde)^2;
gradnormHist = [gradnormHist sqrt(grad_norm)];
fprintf('Iter :%d, norm: %f, gradnorm: %f , regulator: %f\n', iter, norm_iter, grad_norm,regulator);



% normHist(end)
% cost_function_value
% if normHist(end) > cost_function_value || abs((normHist(end) - cost_function_value)/cost_function_value) < 0.2
% if normHist(end) > cost_function_value
    for l = 0:NumOfPoses-2
       AbsolutePoses(:,:,usedPoses(l+2)) = AbsolutePoses_tmp(:,:,usedPoses(l+2));
    end
    FeatureBag(1:3, featureAvailableList+1) = FeatureBag_tmp(1:3, featureAvailableList+1);
%     
%     regulator = regulator * 0.9;
%      
% else
%        regulator = regulator * 1.5;
%        normHist(end) = [];
%        gradnormHist(end) = [];
%        continue;
% end





% grad_norm = norm(f_tilde)^2 + norm(p_tilde)^2;
% gradnormHist = [gradnormHist sqrt(grad_norm)];

% 6. Visualization
% figure(3); clf; hold on;
% Draw3DWorld(AbsolutePoses, FeatureBag, [eye(3) zeros(3,1)], usedPoses);
% grid on;



%     fprintf('Iter :%d, norm: %f, gradnorm: %f , regulator: %f\n', iter, norm_iter, grad_norm,regulator);
    iter = iter + 1;
end
