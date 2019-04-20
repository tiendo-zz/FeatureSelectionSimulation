clear all; close all;clc
% codepath = '/usr/local/google/home/feiwu/Fei/PTAM';

addpath('../utilities/robotics3D/');
addpath('../utilities/PoseGraphMatrix/');
addpath('../utilities/P3P/');
addpath('./pixel/');
addpath('../check_data/');
addpath('../InverseDepth/Triangulation/');
addpath('../InverseDepth/Jacobian/pixel/');
addpath('../InverseDepth/Jacobian/homogeneous/');
addpath('../RobustCostFunction/');
addpath('../TwoViewReconstruction/');
NumOfFeatures = 1000;
NumOfPoses = 50;


CameraMotion = 'StraightToPlane';
%% intrinsic matrix
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
    
%% generate 3d environment
if strcmpi(CameraMotion,'ball')
% [FeatureBag_true_InvDepth, FeatureBag_true_xyz, featureExtracted_true, AbsolutePoses_true, PoseGraphMatrix] = GenerateGroundTruth_ball(NumOfFeatures, NumOfPoses);
[FeatureBag_true_InvDepth, FeatureBag_true_xyz, featureExtracted_true, AbsolutePoses_true, PoseGraphMatrix] = GenerateGroundTruth_ball_2(pi,NumOfFeatures, NumOfPoses,K_true);
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
featureExtracted_true_pixel = cell(size(featureExtracted_true));
for i = 1:length(featureExtracted_true_pixel)
    featureExtracted_true_pixel{i} = K_true*featureExtracted_true{i};
end

%% add noise 
% MaxIter = 3;
% for iter = 0:MaxIter

sigma_homo = 1e-3; % noise for measurement
sigma_q = 0.01; % noise for current camera orientation
sigma_p = 0.1; % 0.2, noise for current camera position
sigma_z = 0.1;

sigma_fc = 30;
sigma_cc = 10;

featureExtracted = cell(NumOfPoses, 1);
f_tilde = cell(NumOfPoses, 1);
for k = 1:NumOfPoses
   f_tilde{k} = sigma_homo*randn(2,size(featureExtracted_true{k},2));
   featureExtracted{k} = [featureExtracted_true{k}(1:2,:) + f_tilde{k};featureExtracted_true{k}(3,:)];
end
AbsolutePoses = zeros(3,4,NumOfPoses);


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
featureExtracted_pixel = cell(size(featureExtracted));
z_tilde = cell(NumOfPoses, 1);
for k = 1:length(featureExtracted_pixel)
    z_tilde{k} = sigma_z*randn(2,size(featureExtracted_true_pixel{k},2));
    featureExtracted_pixel{k} = [featureExtracted_true_pixel{k}(1:2,:) + z_tilde{k};featureExtracted_true_pixel{k}(3,:)];
%    featureExtracted_pixel{k} = CameraParams.K*featureExtracted_true{k};
end

% pertube second pose
dq = [0.5*randn(3,1)*sigma_q;1];
dq = dq / norm(dq);   
C2_R_C1 = quat2rot(quat_mul(dq, rot2quat(AbsolutePoses_true(:,1:3,2))));
C2_t_C1 = -C2_R_C1*(-AbsolutePoses_true(:,1:3,2)'*AbsolutePoses_true(:,4,2) + sigma_p*randn(3,1));
AbsolutePoses(:,:,2) = [C2_R_C1 C2_t_C1];



% save('FeatureBag_true_InvDepth.mat','FeatureBag_true_InvDepth');
% save('FeatureBag_true_xyz.mat','FeatureBag_true_xyz');
% save('featureExtracted_true.mat','featureExtracted_true');
% save('AbsolutePoses_true.mat','AbsolutePoses_true');
% save('PoseGraphMatrix.mat','PoseGraphMatrix');
% save('featureExtracted.mat','featureExtracted');
% save('AbsolutePoses.mat','AbsolutePoses');
% load('FeatureBag_true_InvDepth.mat','FeatureBag_true_InvDepth');
% load('FeatureBag_true_xyz.mat','FeatureBag_true_xyz');
% load('featureExtracted_true.mat','featureExtracted_true');
% load('AbsolutePoses_true.mat','AbsolutePoses_true');
% load('PoseGraphMatrix.mat','PoseGraphMatrix');
% load('featureExtracted.mat','featureExtracted');
% load('AbsolutePoses.mat','AbsolutePoses');

% CameraParams.fc = [1.1 1.1];
% CameraParams.cc = [0.01 0.01];
% CameraParams.K = eye(3);
% CameraParams.Kinv = eye(3);
% [AbsolutePoses_Cartesian,FeaturesBag_Cartesian,~] = simulation_full_BA_Cartesian(PoseGraphMatrix,AbsolutePoses,featureExtracted,NumOfPoses);
%% add outlier
noise_level = 50;
outlier_num = 1*NumOfPoses;
[featureExtracted_outlier,PoseGraphMatrix_outlier] = CreateOutlier(PoseGraphMatrix,featureExtracted_pixel,outlier_num,noise_level);

[AbsolutePoses_InvDep,FeaturesBag_InvDep,~,CameraParams_Invdep] = simulation_full_BA_InvDep_pixel(PoseGraphMatrix,AbsolutePoses,featureExtracted_outlier,NumOfPoses,CameraParams);

error_intrinsic = zeros(4,2);
error_intrinsic(1,1) = abs(CameraParams_Invdep.fc(1) - K_true(1,1));
error_intrinsic(2,1) = abs(CameraParams_Invdep.fc(2) - K_true(2,2));
error_intrinsic(3,1) = abs(CameraParams_Invdep.cc(1) - K_true(1,3));
error_intrinsic(4,1) = abs(CameraParams_Invdep.cc(2) - K_true(2,3));
error_intrinsic(1,2) = abs(CameraParams.fc(1) - K_true(1,1));
error_intrinsic(2,2) = abs(CameraParams.fc(2) - K_true(2,2));
error_intrinsic(3,2) = abs(CameraParams.cc(1) - K_true(1,3));
error_intrinsic(4,2) = abs(CameraParams.cc(2) - K_true(2,3));

figure;
plot(error_intrinsic(:,1),'r*');hold on
plot(error_intrinsic(:,2),'b*');
legend('output','input');
title('error of intrinsic matrix');


% transform FeaturesBag of Inverse Depth param into Cartesian param
FeaturesBag_InvDep_xyz = zeros(size(FeaturesBag_InvDep));
for k = 1:NumOfFeatures
    if FeaturesBag_InvDep(4,k)~=0
       FeaturesBag_InvDep_xyz(1:3,k) =  1/FeaturesBag_InvDep(3,k)*[cos(FeaturesBag_InvDep(2,k))*cos(FeaturesBag_InvDep(1,k));cos(FeaturesBag_InvDep(2,k))*sin(FeaturesBag_InvDep(1,k));sin(FeaturesBag_InvDep(2,k))];
       FeaturesBag_InvDep_xyz(4,k) = FeaturesBag_InvDep(4,k);
    end
end

% reprojection 
for i = 1:NumOfPoses
        subplot(NumOfPoses,1,i); hold on
        plot(featureExtracted_true_pixel{i}(1,:),featureExtracted_true_pixel{i}(2,:),'r*');hold on
%         axis([0 CameraParams.cc(1)*2 0 CameraParams.cc(2)*2]);hold on
        triangulate_feat = find(FeaturesBag_InvDep(4,:) & PoseGraphMatrix(:,i)');
        repro = zeros(3,length(triangulate_feat));
        for k = triangulate_feat
            m_theta_phi = [cos(FeaturesBag_InvDep(2,k))*cos(FeaturesBag_InvDep(1,k));cos(FeaturesBag_InvDep(2,k))*sin(FeaturesBag_InvDep(1,k));sin(FeaturesBag_InvDep(2,k))];
            Ck_P_W = AbsolutePoses_InvDep(:,:,i)*[InversePose(AbsolutePoses_InvDep(:,:,FeaturesBag_InvDep(4,k)));0 0 0 1]*[m_theta_phi;FeaturesBag_InvDep(3,k)];
%             Ck_P_W = AbsolutePoses(:,:,i)*[InversePose(AbsolutePoses(:,:,FeaturesBag(4,k)));0 0 0 1]*[m_theta_phi;1];
            Ck_P_W = Ck_P_W/Ck_P_W(3);
            repro(:,k) = CameraParams_Invdep.K*Ck_P_W;
            plot(repro(1,k),repro(2,k),'bo');hold on
        end        
end

% scale 
% [AbsolutePoses_Cartesian,FeaturesBag_Cartesian,scale_Cartesian] = Scale_Position(AbsolutePoses_true,AbsolutePoses_Cartesian,FeaturesBag_Cartesian,FeatureBag_true_xyz);
[AbsolutePoses_InvDep,FeaturesBag_InvDep_xyz,scale_InvDep] = Scale_Position(AbsolutePoses_true,AbsolutePoses_InvDep,FeaturesBag_InvDep_xyz,FeatureBag_true_xyz);

%% Compute error
% FeaturesBag
error_f = zeros(NumOfFeatures,2);
for k = 1:NumOfFeatures
%     if FeaturesBag_Cartesian(4,k)~=0
%        error_f(k,1) = norm(FeaturesBag_Cartesian(1:3,k) - FeatureBag_true_xyz(1:3,k));
%     end
    if FeaturesBag_InvDep(4,k)~=0
       error_f(k,2) = norm(FeaturesBag_InvDep_xyz(1:3,k) - FeatureBag_true_xyz(1:3,k));
    end
end
figure;
plot(error_f(:,1));hold on
plot(error_f(:,2));
legend('Cartesian','Inverse Depth');
title(['feature error, sigma z : ',num2str(sigma_homo)]);
% saveas(gcf,['result/',CameraMotion,'/Feature_2viewLS_',CameraMotion,'_error with sigma_z : ',num2str(sigma_z),'.fig']);
% saveas(gcf,['result/',CameraMotion,'/Feature_2viewLS_',CameraMotion,'_error with sigma_z : ',num2str(sigma_z),'.jpg']);

% camera pose
error_r = zeros(NumOfPoses,3);
error_t = zeros(NumOfPoses,3);
for k = 1:NumOfPoses   
    error_r(k,1) = acos((trace(AbsolutePoses_Cartesian(:,1:3,k)*AbsolutePoses_true(:,1:3,k)')-1)/2);
    error_t(k,1) = sum((AbsolutePoses_Cartesian(:,4,k)- AbsolutePoses_true(:,4,k)).^2)^0.5;
    error_r(k,2) = acos((trace(AbsolutePoses_InvDep(:,1:3,k)*AbsolutePoses_true(:,1:3,k)')-1)/2);
    error_t(k,2) = sum((AbsolutePoses_InvDep(:,4,k)- AbsolutePoses_true(:,4,k)).^2)^0.5;
end

figure;
subplot(2,1,1);
plot(error_r(:,1));hold on
plot(error_r(:,2));
legend('Cartesian','Inverse Depth');
title('rotation');
subplot(2,1,2);
plot(error_t(:,1));hold on
plot(error_t(:,2));
legend('Cartesian','Inverse Depth');
title(['position, sigma z : ',num2str(sigma_homo)]);
% saveas(gcf,['result/',CameraMotion,'/CameraPose_2viewLS_',CameraMotion,'_error with sigma_z : ',num2str(sigma_z),'.fig']);
% saveas(gcf,['result/',CameraMotion,'/CameraPose_2viewLS_',CameraMotion,'_error with sigma_z : ',num2str(sigma_z),'.jpg']);