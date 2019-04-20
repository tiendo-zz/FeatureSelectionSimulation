function [FeatureBag, AbsolutePoses, CameraParams,PoseGraphMatrix,normHist, gradnormHist,EraseId,error] = BA_InverseDepth_2viewBLS_Robust_pixel_simulation(FeatureBag, AbsolutePoses, PoseGraphMatrix, featureExtracted, usedPoses, CameraParams,sigma)


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
maxIter = 50;
normHist = zeros(maxIter,1);
gradnormHist = zeros(maxIter,1);
iter = 1;
thresh = 30;
sigma = 2*thresh;

% 1. Build Jacobian Jp Jf corresponds to each feature
featureAvailableList_all = find(FeatureBag(4,:));

featureAvailableList = find(FeatureBag(4,:));
featureAvailableList = featureAvailableList-1;


ErasePose = cell(featureAvailableList_all(end),1);
for e = 1:length(ErasePose)
    ErasePose{e} = [];
end


while iter < maxIter && norm_iter > 1 && grad_norm > 1e-4
    
AbsolutePoses_tmp = AbsolutePoses;
FeatureBag_tmp = FeatureBag;
   
NumOfFeatures = length(featureAvailableList);    
poses_see_featk_all = length(find(PoseGraphMatrix_r(featureAvailableList+1, :)));


feat_cam = zeros(poses_see_featk_all,2);

Jpose = zeros(2*poses_see_featk_all, 6*(NumOfPoses-1)+5);
Jfc_cc = zeros(2*poses_see_featk_all, 4);
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
      Ck_z_f = featureExtracted{usedPoses(poses_see_featk(l+1))}(1:2,PoseGraphMatrix_r(featk+1,poses_see_featk(l+1)));                  
      
      
      % Pose Transformation
      Ck_T_W = AbsolutePoses(:,:,usedPoses(poses_see_featk(l+1)));
      Cr_T_W = AbsolutePoses(:,:,FeatureBag(4, featk+1));          
      
      Ck_R_W = Ck_T_W(:,1:3);
      W_R_Cr = Cr_T_W(:,1:3)';
      W_p_Cr = -Cr_T_W(:,1:3)'*Cr_T_W(:,4);
      W_p_Ck = -Ck_T_W(:,1:3)'*Ck_T_W(:,4);
      
      
      m_theta_phi = [cos(Cr_p_f(2))*cos(Cr_p_f(1));cos(Cr_p_f(2))*sin(Cr_p_f(1));sin(Cr_p_f(2))];
  
      ref_idx_usePoses = find(FeatureBag(4, featk+1) == usedPoses);
      
      Ck_p_f = Ck_R_W*W_R_Cr*m_theta_phi + Cr_p_f(3)*Ck_R_W*(W_p_Cr-W_p_Ck);
      H_pixel = [CameraParams.fc(1) 0;0 CameraParams.fc(2)];
      Ck_z_f_hat = CameraParams.K*(Ck_p_f / Ck_p_f(3));
      
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
         residual_tmp =  Ck_z_f - Ck_z_f_hat(1:2,:);  
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
          [Jq,~,~] = dprojection_dq(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);          
          Jqr = dprojection_dqr(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
          Jro = dprojection_dro(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
      end
      Jtheta = dprojection_dtheta(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
      Jphi = dprojection_dphi(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f);
      
      [Jfc,Jcc,residual_tmp,flag] = dprojection_dfc_dcc_simulation(Cr_p_f, Ck_R_W, W_R_Cr,W_p_Cr, W_p_Ck, Ck_z_f,CameraParams,thresh);
      if isempty(Jfc)
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
      
      residual_r = Robust_CostFunction(residual_tmp,sigma);
      Jr = Robust_Cost_Jacobian(residual_tmp,sigma); 
      
      if usedPoses(poses_see_featk(l+1)) == 2
         Jpose(2*meascount+1:2*meascount+2, 6*(poses_see_featk(l+1)-1)+1:6*(poses_see_featk(l+1)-1)+5) = Jr*H_pixel*[Jp Jq];
      else
          if usedPoses(poses_see_featk(l+1)) ~= 1      
              Jpose(2*meascount+1:2*meascount+2, 6*(poses_see_featk(l+1)-2)+5+1:6*(poses_see_featk(l+1)-2)+5+6) = Jr*H_pixel*[Jp Jq];
          end
      end
      
      
      if ref_idx_usePoses ~= 1      
         if usedPoses(ref_idx_usePoses) == 2
            Jpose(2*meascount+1:2*meascount+2, 6*(ref_idx_usePoses-1)+1:6*(ref_idx_usePoses-1)+5) = Jr*H_pixel*[Jpr Jqr];
         else
             if usedPoses(ref_idx_usePoses) ~= 1      
                Jpose(2*meascount+1:2*meascount+2, 6*(ref_idx_usePoses-2)+5+1:6*(ref_idx_usePoses-2)+5+6) = Jr*H_pixel*[Jpr Jqr];
             end
         end
      end
      
      residual(2*meascount+1:2*meascount+2) = residual_r;
      Jfeat(2*meascount+1:2*meascount+2, 3*featcount+1:3*featcount+3) = Jr*H_pixel*[Jtheta Jphi Jro];              
      
      Jfc_cc(2*meascount+1:2*meascount+2, :) = Jr*[Jfc Jcc];
      
      feat_cam(meascount+1,:) = [featk+1 poses_see_featk(l+1)];
      meascount = meascount + 1;
      count_measure_single = count_measure_single + 1;
   end   
   
   if ((count_measure_single < 3) && (NumOfPoses > 5)) || count_measure_single == 1
       outlier_num = outlier_num + 1;
       EraseId(outlier_num) = featk;
       for i = meascount : meascount - count_measure_single
           Jpose(2*i+1:2*i+2, :) = 0;
           Jfc_cc(2*i+1:2*i+2, :) = 0;
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
EraseId(outlier_num+1:end) = []; 
Jpose(:,1:6) = [];


% clear the all zero lines
Jfeat(:,3*featcount+1:end) = [];
Jfeat(2*meascount+1:end,:) = [];
Jpose(2*meascount+1:end,:) = [];
Jfc_cc(2*meascount+1:end,:) = [];
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

% Jpose_sp = sparse([Jpose Jfc_cc]); %Jpose_sp contains both Jpose and Jacobian of extrinsic matrix
% Jfeat_sp = sparse(Jfeat);
% A11 = Jpose_sp'*Jpose_sp;
% A12 = Jpose_sp'*Jfeat_sp;
% A22_inv_sp = sparse(A22_inv);
% res_ = [Jpose_sp'*residual;Jfeat_sp'*residual];
Jpose_sp = sparse(Jpose); 
Jfc_cc_sp = sparse(Jfc_cc);
Jfeat_sp = sparse(Jfeat);
A11 = [Jpose_sp'*Jpose_sp Jpose_sp'*Jfc_cc_sp;Jfc_cc_sp'*Jpose_sp Jfc_cc_sp'*Jfc_cc_sp];
A12 = [Jpose_sp'*Jfeat_sp;Jfc_cc_sp'*Jfeat_sp];
A22_inv_sp = sparse(A22_inv);
res_ = [Jpose_sp'*residual;Jfc_cc_sp'*residual;Jfeat_sp'*residual];

% solve:
p_tilde = (A11 +regulator*speye(size(A11))- A12*A22_inv_sp*A12') \ (res_(1:6*(NumOfPoses-2)+5+4) - A12*A22_inv_sp*res_(1+6*(NumOfPoses-2)+5+4:end));

% do not update first 3 camera poses (if more than 3 poses)
for l = 0:NumOfPoses-2
   if l == 0
       dq = [0.5*p_tilde(6*l+3:6*l+5);1];
       dq = dq / norm(dq);
       c1_t_c2 = -AbsolutePoses_tmp(:,1:3,usedPoses(l+2))'*AbsolutePoses_tmp(:,4,usedPoses(l+2));
       AbsolutePoses_tmp(:,1:3,usedPoses(l+2)) = quat2rot( quat_mul(dq, rot2quat(AbsolutePoses_tmp(:,1:3,usedPoses(l+2)))) );      
       c1_t_c2 = c1_t_c2 - p_tilde(6*l+1)*t_p + p_tilde(6*l+2)*t_pp;
%        c1_t_c2 = c1_t_c2 / norm(c1_t_c2);
       AbsolutePoses_tmp(:,4,usedPoses(l+2)) = -AbsolutePoses_tmp(:,1:3,usedPoses(l+2))*c1_t_c2;
   else 
       dq = [0.5*p_tilde(6*(l-1)+5+4:6*(l-1)+5+6);1];
       dq = dq / norm(dq);
       W_p_Ck = -AbsolutePoses_tmp(:,1:3,usedPoses(l+2))'*AbsolutePoses_tmp(:,4,usedPoses(l+2)) + p_tilde(6*(l-1)+5+1:6*(l-1)+5+3);
       AbsolutePoses_tmp(:,1:3,usedPoses(l+2)) = quat2rot( quat_mul(dq, rot2quat(AbsolutePoses_tmp(:,1:3,usedPoses(l+2)))) );       
       AbsolutePoses_tmp(:,4,usedPoses(l+2)) = -AbsolutePoses_tmp(:,1:3,usedPoses(l+2))*W_p_Ck;
   end
end


% update intrinsic matrix
CameraParams.fc = CameraParams.fc + [p_tilde(end-3) p_tilde(end-2)];
CameraParams.cc = CameraParams.cc + [p_tilde(end-1) p_tilde(end)];
fx = CameraParams.fc(1); fy = CameraParams.fc(2);
cx = CameraParams.cc(1); cy = CameraParams.cc(2);
CameraParams.K = [fx 0 cx; ...
                  0 fy cy; ...
                  0 0 1];
CameraParams.Kinv = [1/fx 0 -cx/fx; ...
                     0 1/fy -cy/fy; ...
                     0 0 1]; 
% 4. Efficient solve for each feature f
bf_p = res_(1+6*(NumOfPoses-2)+5+4:end) - A12'*p_tilde;
f_tilde = A22_inv_sp * bf_p;
[~,erase_b] = ismember(EraseId,featureAvailableList);
FeatureBag(:, featureAvailableList(erase_b)+1) = 0; % erase outliers
featureAvailableList(erase_b) = [];


FeatureBag_tmp(1:3, featureAvailableList+1) = FeatureBag_tmp(1:3, featureAvailableList+1) + reshape(f_tilde, [3, length(featureAvailableList)]);

%%after updating the temp result, compute the cost function one more time,
%%if it's smaller than the previous, take it and decrease regulator, else,
%%discard it and increase twice the regulator



grad_norm = sqrt(norm(f_tilde)^2 + norm(p_tilde)^2);
gradnormHist(iter) = grad_norm;
fprintf('Iter :%d, norm: %f, gradnorm: %f , regulator: %f, thresh: %f\n', iter, norm_iter, grad_norm,regulator,thresh);

if iter >= 2
    if (norm_iter > normHist(iter-1)) && (grad_norm > gradnormHist(iter-1))
        thresh = thresh * 0.8;
        sigma = thresh*2;
    elseif (norm_iter < normHist(iter-1)) && (grad_norm < gradnormHist(iter-1))
        thresh = thresh * 1.25;
        sigma = thresh*2;
    end
end

for l = 0:NumOfPoses-2
    AbsolutePoses(:,:,usedPoses(l+2)) = AbsolutePoses_tmp(:,:,usedPoses(l+2));
end
FeatureBag(1:3, featureAvailableList+1) = FeatureBag_tmp(1:3, featureAvailableList+1);

iter = iter + 1;
    
    
end

normHist(iter:end) = [];
gradnormHist(iter:end) = [];
