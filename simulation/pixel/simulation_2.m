clear all; close all;clc
codepath = '/usr/local/google/home/feiwu/Fei/PTAM_code';

addpath([codepath,'/simulation/']);
addpath([codepath,'/utilities/robotics3D/']);
addpath([codepath,'/utilities/PoseGraphMatrix/']);
addpath([codepath,'/utilities/P3P/']);
addpath([codepath,'/check_data/']);
addpath([codepath,'/InverseDepth/Triangulation/']);
addpath([codepath,'/InverseDepth/Jacobian/pixel/']);
addpath([codepath,'/RobustCostFunction/']);
addpath([codepath,'/TwoViewReconstruction/']);
addpath([codepath,'/EPnP_matlab/EPnP/']);
addpath([codepath,'/EPnP_matlab/error/']);
NumOfFeatures = 100;
NumOfPoses = 5;


CameraMotion = 'ball';

if strcmpi(CameraMotion,'ball')
[FeatureBag_true_InvDepth, FeatureBag_true_xyz, featureExtracted_true, AbsolutePoses_true, PoseGraphMatrix] = GenerateGroundTruth_ball(NumOfFeatures, NumOfPoses);
elseif strcmpi(CameraMotion,'circle')
[FeatureBag_true_InvDepth, FeatureBag_true_xyz, featureExtracted_true, AbsolutePoses_true, PoseGraphMatrix] = GenerateGroundTruth_Circle(NumOfFeatures, NumOfPoses);
elseif strcmpi(CameraMotion,'StraightSideWay')
[FeatureBag_true_InvDepth, FeatureBag_true_xyz, featureExtracted_true, AbsolutePoses_true, PoseGraphMatrix] = GenerateGroundTruth_StraightSideWay(NumOfFeatures, NumOfPoses);
elseif strcmpi(CameraMotion,'StraightToPlane')
[FeatureBag_true_InvDepth, FeatureBag_true_xyz, featureExtracted_true, AbsolutePoses_true, PoseGraphMatrix] = GenerateGroundTruth_StraightToPlane(NumOfFeatures, NumOfPoses);
end
VisualizeMultiPoses(AbsolutePoses_true, FeatureBag_true_xyz, 1:(NumOfPoses-1), NumOfPoses);
% saveas(gcf,['result/',CameraMotion,'/3d_environment_',CameraMotion,'.fig']);

%% convert to pixel
fc = [254.997, 254.933];
cc = [326.69, 230.118];
fx = fc(1); fy = fc(2);
cx = cc(1); cy = cc(2);
K_true = [fx 0 cx; ...
        0 fy cy; ...
        0 0 1];
Kinv_true = [1/fx 0 -cx/fx; ...
        0 1/fy -cy/fy; ...
        0 0 1];
featureExtracted_true_pixel = cell(size(featureExtracted_true));
for i = 1:length(featureExtracted_true_pixel)
    featureExtracted_true_pixel{i} = K_true*featureExtracted_true{i};
end

%% add noise 
% MaxIter = 3;
% for iter = 0:MaxIter

sigma_homo = 1e-3; % noise for measurement
sigma_q = 0.1; % noise for current camera orientation
sigma_p = 0.1; % 0.2, noise for current camera position
sigma_z = 0.1;



featureExtracted = cell(NumOfPoses, 1);
f_tilde = cell(NumOfPoses, 1);
for k = 1:NumOfPoses
   f_tilde{k} = sigma_homo*randn(2,size(featureExtracted_true{k},2));
   featureExtracted{k} = [featureExtracted_true{k}(1:2,:) + f_tilde{k};featureExtracted_true{k}(3,:)];
end

featureExtracted_pixel = cell(size(featureExtracted));
z_tilde = cell(NumOfPoses, 1);
for k = 1:length(featureExtracted_pixel)
    z_tilde{k} = sigma_z*randn(2,size(featureExtracted_true_pixel{k},2));
    featureExtracted_pixel{k} = [featureExtracted_true_pixel{k}(1:2,:) + z_tilde{k};featureExtracted_true_pixel{k}(3,:)];
%    featureExtracted_pixel{k} = CameraParams.K*featureExtracted_true{k};
end

% pertube poses
AbsolutePoses = AbsolutePoses_true;
for k = 3:NumOfPoses
    dq = [0.5*randn(3,1)*sigma_q;1];
    dq = dq / norm(dq);   
    Ck_R_C1 = quat2rot(quat_mul(dq, rot2quat(AbsolutePoses_true(:,1:3,k))));
    Ck_t_C1 = -Ck_R_C1*(-AbsolutePoses_true(:,1:3,k)'*AbsolutePoses_true(:,4,k) + sigma_p*randn(3,1));
    AbsolutePoses(:,:,k) = [Ck_R_C1 Ck_t_C1];
end


iter = 1;
maxIter = 5;
error_t = zeros(NumOfPoses,maxIter);

while iter < maxIter

sigma_fc = iter*10;
sigma_cc = 1;

fc = fc+rand(1,2)*sigma_fc;
cc = cc+rand(1,2)*sigma_cc;
fx = fc(1); fy = fc(2);
cx = cc(1); cy = cc(2);
CameraParams.K = [fx 0 cx; ...
        0 fy cy; ...
        0 0 1];
CameraParams.Kinv = [1/fx 0 -cx/fx; ...
        0 1/fy -cy/fy; ...
        0 0 1];
CameraParams.fc = fc; CameraParams.cc = cc;


usedPoses  = [1;2];
AbsolutePoses_InvDep(:,:,1:2) = AbsolutePoses_true(:,:,1:2);
for k = 3:NumOfPoses
    InlierIdx = find(PoseGraphMatrix(:,k));
    
    [Cn_R_Cr, Cn_t_Cr, normHist] = simulation_PnP_NL_InvDep_pixel(AbsolutePoses(:,1:3,k), AbsolutePoses(:,4,k),AbsolutePoses, FeatureBag_true_InvDepth(:, InlierIdx), featureExtracted_pixel{k}(1:2,:),CameraParams);                  
    AbsolutePoses(:,:,k) = [Cn_R_Cr, Cn_t_Cr];
      

%     usedPoses = [usedPoses;k];
%     [FeatureBag_InvDep, AbsolutePoses_InvDep, CameraParams_InvDep,~, ~,~] = BA_Intrinsic_simulation(FeatureBag_true_InvDepth, AbsolutePoses, PoseGraphMatrix, featureExtracted_pixel, usedPoses, CameraParams);
end

error_intrinsic = zeros(4,2);
error_intrinsic(1,1) = abs(CameraParams_InvDep.fc(1) - K_true(1,1));
error_intrinsic(2,1) = abs(CameraParams_InvDep.fc(2) - K_true(2,2));
error_intrinsic(3,1) = abs(CameraParams_InvDep.cc(1) - K_true(1,3));
error_intrinsic(4,1) = abs(CameraParams_InvDep.cc(2) - K_true(2,3));
error_intrinsic(1,2) = abs(CameraParams.fc(1) - K_true(1,1));
error_intrinsic(2,2) = abs(CameraParams.fc(2) - K_true(2,2));
error_intrinsic(3,2) = abs(CameraParams.cc(1) - K_true(1,3));
error_intrinsic(4,2) = abs(CameraParams.cc(2) - K_true(2,3));
figure;
plot(error_intrinsic(:,1),'r*');hold on
plot(error_intrinsic(:,2),'b*');
legend('output','input');
title('error of intrinsic matrix');

% % camera pose
% figure;
% error_r = zeros(NumOfPoses,3);
% error_t = zeros(NumOfPoses,3);
% for k = 1:NumOfPoses   
%     error_r(k,2) = acos((trace(AbsolutePoses_InvDep(:,1:3,k)*AbsolutePoses_true(:,1:3,k)')-1)/2);
%     error_t(k,2) = sum((AbsolutePoses_InvDep(:,4,k)- AbsolutePoses_true(:,4,k)).^2)^0.5;
% end
% plot(error_r(:,2));hold on
% plot(error_t(:,2));
% legend('rotation','position');

for k = 1:NumOfPoses   
    error_t(k,iter) = sum((AbsolutePoses_InvDep(:,4,k)- AbsolutePoses_true(:,4,k)).^2)^0.5;
end
iter = iter + 1;
end

figure;
for iter = 1:maxIter
    plot(error_t(:,iter)); hold on
end
legend(['sigma fc = ',num2str(10)],['sigma fc = ',num2str(20)],['sigma fc = ',num2str(30)],['sigma fc = ',num2str(40)],['sigma fc = ',num2str(50)]);