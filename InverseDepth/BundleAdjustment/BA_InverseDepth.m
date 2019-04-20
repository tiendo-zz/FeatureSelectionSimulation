function [FeatureBag, AbsolutePoses, normHist, gradnormHist,EraseId,error] = BA_InverseDepth(FeatureBag, AbsolutePoses, PoseGraphMatrix, featureExtracted, usedPoses, CameraParams)

PoseGraphMatrix_r = PoseGraphMatrix(:, usedPoses);
NumOfPoses = length(usedPoses);


%% Bundle Adjustment
regulator = 0.1;
norm_iter = 100;
grad_norm = 100;
maxIter = 50;
normHist = zeros(maxIter,1);
gradnormHist = zeros(maxIter,1);
iter = 1;

    
% 1. Build Jacobian Jp Jf corresponds to each feature
featureAvailableList_all = find(FeatureBag(4,:));

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

Jpose = zeros(2*poses_see_featk_all, 6*NumOfPoses);
Jfeat = zeros(2*poses_see_featk_all, 3*NumOfFeatures);
error = zeros(poses_see_featk_all,2);
residual = zeros(2*poses_see_featk_all, 1);
meascount = 0;
last_meas_count = 0;

featcount = 0;
% Hessian build up
A22_inv = zeros(3*NumOfFeatures);


EraseId = zeros(length(featureAvailableList),1);
outlier_num = 0;
for featk = featureAvailableList %0:size(PoseGraphMatrix,1)-1
   
    
   poses_see_featk = find(PoseGraphMatrix_r(featk+1, :));
   if ~isempty(ErasePose{featk+1})
       [~,b] = ismember(ErasePose{featk+1},poses_see_featk);
       poses_see_featk(b) = [];
   end
   
   if (NumOfPoses > 5) && (length(poses_see_featk) < 3)
       outlier_num = outlier_num + 1;
       EraseId(outlier_num) = featk;
       continue;
   end
   
   count_measure_single = 0;
   for l = 0:length(poses_see_featk)-1
       
      %  Pose Jacobian      
      Cr_p_f = FeatureBag(1:3 ,featk+1);      
      Ck_z_f = CameraParams.Kinv(1:2,:)*[featureExtracted{usedPoses(poses_see_featk(l+1))}(1:2,PoseGraphMatrix_r(featk+1,poses_see_featk(l+1)));1];                  
      
      
      % Pose Transformation
      Ck_T_W = AbsolutePoses(:,:,usedPoses(poses_see_featk(l+1)));
      Cr_T_W = AbsolutePoses(:,:,FeatureBag(4, featk+1));                
      
      Ck_R_W = Ck_T_W(:,1:3);
      W_R_Cr = Cr_T_W(:,1:3)';
      W_p_Cr = -Cr_T_W(:,1:3)'*Cr_T_W(:,4);
      W_p_Ck = -Ck_T_W(:,1:3)'*Ck_T_W(:,4);
      
      
      m_theta_phi = [cos(Cr_p_f(2))*cos(Cr_p_f(1));cos(Cr_p_f(2))*sin(Cr_p_f(1));sin(Cr_p_f(2))];
      Ck_z_f_hat = m_theta_phi(1:2) / m_theta_phi(3);
      
      if usedPoses(poses_see_featk(l+1))==FeatureBag(4, featk+1)
         Jq = zeros(2,3);
         Jp = zeros(2,3);
         Jqr = zeros(2,3);
         Jpr = zeros(2,3);
         Jro = zeros(2,1);
         residual_tmp = Ck_z_f - Ck_z_f_hat;
      else
          [Jq,residual_tmp,flag] = dprojection_dq(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
          Jp = dprojection_dp(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
          Jqr = dprojection_dqr(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
          Jpr = dprojection_dpr(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
          Jro = dprojection_dro(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
      end
      Jtheta = dprojection_dtheta(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
      Jphi = dprojection_dphi(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
      
      if isempty(Jq)
         % if manhalanobis dist is twice greater than threshold, delete
         % this measurement (flag = 2); otherwise, record this camera pose
         % idx (flag = 1).
         if flag == 1
            ErasePose{featk+1} = [ErasePose{featk+1};poses_see_featk(l+1)];
         elseif flag == 2
            PoseGraphMatrix(featk+1,usedPoses(poses_see_featk(l+1))) = 0;
            PoseGraphMatrix_r(featk+1,poses_see_featk(l+1)) = 0;
         end
         continue;
      end
      
           

      Jpose(2*meascount+1:2*meascount+2, 6*(poses_see_featk(l+1)-1)+1:6*(poses_see_featk(l+1)-1)+6) = [Jp Jq];      
      ref_idx_usePoses = find(FeatureBag(4, featk+1) == usedPoses);
      Jpose(2*meascount+1:2*meascount+2, 6*(ref_idx_usePoses-1)+1:6*(ref_idx_usePoses-1)+6) = [Jpr Jqr];      
      
      residual(2*meascount+1:2*meascount+2) = residual_tmp;
      
      Jfeat(2*meascount+1:2*meascount+2, 3*featcount+1:3*featcount+3) = [Jtheta Jphi Jro];              
      
      feat_cam(meascount+1,:) = [featk+1 poses_see_featk(l+1)];
      meascount = meascount + 1;
      count_measure_single = count_measure_single + 1;
   end   
   
   if (count_measure_single < 3) && (NumOfPoses > 5)
       outlier_num = outlier_num + 1;
       EraseId(outlier_num) = featk;
       for i = meascount : meascount - count_measure_single
           Jpose(2*i+1:2*i+2, :) = 0;
           residual(2*i+1:2*i+2) = 0;
           Jfeat(2*i+1:2*i+2, 3*featcount+1:3*featcount+3) = 0;              
           feat_cam(i+1,:) = 0;
       end
       meascount = meascount - count_measure_single;
       continue;
   end
   
   Jfeatk = Jfeat(last_meas_count+1:2*meascount, 3*featcount+1:3*featcount+3);
   A22_inv(3*featcount+1:3*featcount+3, 3*featcount+1:3*featcount+3) = (Jfeatk'*Jfeatk + regulator*eye(3)) \ eye(3);
      
   featcount = featcount + 1;
   last_meas_count = 2*meascount;
end

Jpose(:,1:12) = [];


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
normHist(iter) = norm_iter;


Jpose_sp = sparse(Jpose);
Jfeat_sp = sparse(Jfeat);
A11 = Jpose_sp'*Jpose_sp;
A12 = Jpose_sp'*Jfeat_sp;
A22_inv_sp = sparse(A22_inv);
res_ = [Jpose_sp'*residual;Jfeat_sp'*residual];

% solve:
p_tilde = (A11 +regulator*speye(size(A11))- A12*A22_inv_sp*A12') \ (res_(1:6*(NumOfPoses-2)) - A12*A22_inv_sp*res_(1+6*(NumOfPoses-2):end));

% do not update first 2 camera poses
   for l = 0:NumOfPoses-3
       dq = [0.5*p_tilde(6*l+4:6*l+6);1];
       dq = dq / norm(dq);
       W_p_Ck = -AbsolutePoses_tmp(:,1:3,usedPoses(l+3))'*AbsolutePoses_tmp(:,4,usedPoses(l+3)) + p_tilde(6*l+1:6*l+3);
       AbsolutePoses_tmp(:,1:3,usedPoses(l+3)) = quat2rot( quat_mul(dq, rot2quat(AbsolutePoses_tmp(:,1:3,usedPoses(l+3)))) );       
       AbsolutePoses_tmp(:,4,usedPoses(l+3)) = -AbsolutePoses_tmp(:,1:3,usedPoses(l+3))*W_p_Ck;
   end

 
% 4. Efficient solve for each feature f
bf_p = res_(1+6*(NumOfPoses-2):end) - A12'*p_tilde;
f_tilde = A22_inv_sp * bf_p;
[~,erase_b] = ismember(EraseId,featureAvailableList);
% FeatureBag(1:3, featureAvailableList(erase_b)+1) = 0; % erase outliers
featureAvailableList(erase_b) = [];
FeatureBag_tmp(1:3, featureAvailableList+1) = FeatureBag_tmp(1:3, featureAvailableList+1) + reshape(f_tilde, [3, length(featureAvailableList)]);

%%after updating the temp result, compute the cost function one more time,
%%if it's smaller than the previous, take it and decrease regulator, else,
%%discard it and increase twice the regulator



grad_norm = norm(f_tilde)^2 + norm(p_tilde)^2;
gradnormHist(iter) = sqrt(grad_norm);
fprintf('Iter :%d, norm: %f, gradnorm: %f , regulator: %f\n', iter, norm_iter, grad_norm,regulator);



    for l = 0:NumOfPoses-3
       AbsolutePoses(:,:,usedPoses(l+3)) = AbsolutePoses_tmp(:,:,usedPoses(l+3));
    end
    FeatureBag(1:3, featureAvailableList+1) = FeatureBag_tmp(1:3, featureAvailableList+1);

    iter = iter + 1;
end

normHist(iter:end) = [];
gradnormHist(iter:end) = [];

