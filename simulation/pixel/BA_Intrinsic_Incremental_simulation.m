function [FeatureBag, AbsolutePoses, CameraParams,normHist, gradnormHist,error] = BA_Intrinsic_Incremental_simulation(FeatureBag, AbsolutePoses, PoseGraphMatrix, featureExtracted, usedPoses, CameraParams)

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

while iter < maxIter && norm_iter > 1e-2 && grad_norm > 1e-5
    
AbsolutePoses_tmp = AbsolutePoses;
FeatureBag_tmp = FeatureBag;
   
NumOfFeatures = length(featureAvailableList);    
poses_see_featk_all = length(find(PoseGraphMatrix_r(featureAvailableList+1, :)));


feat_cam = zeros(poses_see_featk_all,2);

Jfc_cc = zeros(2*poses_see_featk_all, 4);
Jfx_cx = zeros(2*poses_see_featk_all, 2);
Jfy_cy = zeros(2*poses_see_featk_all, 2);
error = zeros(poses_see_featk_all,2);
residual = zeros(2*poses_see_featk_all, 1);
meascount = 0;
last_meas_count = 0;

featcount = 0;


for featk = featureAvailableList %0:size(PoseGraphMatrix,1)-1
   
    
   poses_see_featk = find(PoseGraphMatrix_r(featk+1, :));
   
   count_measure_single = 0;
   for l = 0:length(poses_see_featk)-1
       
      %  Pose Jacobian      
      Cr_p_f = FeatureBag(1:3 ,featk+1);      
      Ck_z_f = featureExtracted{usedPoses(poses_see_featk(l+1))}(1:2,PoseGraphMatrix_r(featk+1,poses_see_featk(l+1)));                  
      
      
      % Pose Transformation
      Ck_T_W = AbsolutePoses(:,:,usedPoses(poses_see_featk(l+1)));
      Cr_T_W = AbsolutePoses(:,:,FeatureBag(4, featk+1));               
      
      Ck_R_W = Ck_T_W(:,1:3);
      W_R_Cr = Cr_T_W(:,1:3)';
      W_p_Cr = -Cr_T_W(:,1:3)'*Cr_T_W(:,4);
      W_p_Ck = -Ck_T_W(:,1:3)'*Ck_T_W(:,4);
      
      
      m_theta_phi = [cos(Cr_p_f(2))*cos(Cr_p_f(1));cos(Cr_p_f(2))*sin(Cr_p_f(1));sin(Cr_p_f(2))];
%       Ck_z_f_hat = m_theta_phi(1:2) / m_theta_phi(3);
      
      
      Ck_p_f = Ck_R_W*W_R_Cr*m_theta_phi + Cr_p_f(3)*Ck_R_W*(W_p_Cr-W_p_Ck);
    
      homo = Ck_p_f(1:2) / Ck_p_f(3);
      Ck_z_f_hat = CameraParams.K*[homo;1];
    
      % residual
      residual_tmp = Ck_z_f - Ck_z_f_hat(1:2);
      % Jacobian
      Jfc = [homo(1) 0;0 homo(2)];  
      Jcc = eye(2);
      
      
      residual(2*meascount+1:2*meascount+2) = residual_tmp;            
      
      Jfc_cc(2*meascount+1:2*meascount+2, :) = [Jfc Jcc];
      Jfx_cx(2*meascount+1:2*meascount+2, :) = [Jfc(:,1) Jcc(:,1)];
      Jfy_cy(2*meascount+1:2*meascount+2, :) = [Jfc(:,2) Jcc(:,2)];
      feat_cam(meascount+1,:) = [featk+1 poses_see_featk(l+1)];
      meascount = meascount + 1;
      count_measure_single = count_measure_single + 1;
   end   
      
   featcount = featcount + 1;
   last_meas_count = 2*meascount;
end


% clear the all zero lines
Jfc_cc(2*meascount+1:end,:) = [];
Jfx_cx(2*meascount+1:end,:) = [];
Jfy_cy(2*meascount+1:end,:) = [];
residual(2*meascount+1:end,:) = [];

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

A11 = Jfc_cc'*Jfc_cc;
dx = ( A11 + regulator*speye(size(A11))) \ (Jfc_cc'*residual);
% A11 = Jfx_cx'*Jfx_cx;
% dx = ( A11 + regulator*speye(size(A11))) \ (Jfx_cx'*residual);
% A11 = Jfy_cy'*Jfy_cy;
% dx = ( A11 + regulator*speye(size(A11))) \ (Jfy_cy'*residual);


% update intrinsic matrix
CameraParams.fc(1) = CameraParams.fc(1) + dx(1);
CameraParams.fc(2) = CameraParams.fc(2) + dx(2);
CameraParams.cc(1) = CameraParams.cc(1) + dx(3);
CameraParams.cc(2) = CameraParams.cc(2) + dx(4);
% CameraParams.fc(2) = CameraParams.fc(2) + dx(1);
% CameraParams.cc(2) = CameraParams.cc(2) + dx(2);
fx = CameraParams.fc(1); fy = CameraParams.fc(2);
cx = CameraParams.cc(1); cy = CameraParams.cc(2);
CameraParams.K = [fx 0 cx; ...
                  0 fy cy; ...
                  0 0 1];
CameraParams.Kinv = [1/fx 0 -cx/fx; ...
                     0 1/fy -cy/fy; ...
                     0 0 1]; 
                 

grad_norm = norm(dx)^2;
gradnormHist = [gradnormHist sqrt(grad_norm)];
fprintf('Iter :%d, norm: %f, gradnorm: %f , regulator: %f\n', iter, norm_iter, grad_norm,regulator);

    iter = iter + 1;
end
