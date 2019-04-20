clear all; close all;clc
codepath = '/usr/local/google/home/feiwu/Fei/PTAM_code';

addpath([codepath,'/simulation/']);
addpath([codepath,'/utilities/robotics3D/']);
addpath([codepath,'/utilities/PoseGraphMatrix/']);
addpath([codepath,'/utilities/P3P/']);
addpath([codepath,'/check_data/']);
addpath([codepath,'/InverseDepth/Triangulation/']);
addpath([codepath,'/InverseDepth/Jacobian/pixel/']);
addpath([codepath,'/InverseDepth/Jacobian/homogeneous/']);
addpath([codepath,'/RobustCostFunction/']);
addpath([codepath,'/TwoViewReconstruction/']);
addpath([codepath,'/EPnP_matlab/']);
addpath([codepath,'/EPnP_matlab/EPnP/']);
addpath([codepath,'/EPnP_matlab/error/']);
NumOfFeatures = 100;
NumOfPoses = 50;


CameraMotion = 'ball';

%% intriinsic matrix 
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
sigma_q = 0.05; % noise for current camera orientation
sigma_p = 0.1; % 0.2, noise for current camera position
sigma_z = 1;



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
for k = 2:NumOfPoses
    dq = [0.5*randn(3,1)*sigma_q;1];
    dq = dq / norm(dq);   
    Ck_R_C1 = quat2rot(quat_mul(dq, rot2quat(AbsolutePoses_true(:,1:3,k))));
    Ck_t_C1 = -Ck_R_C1*(-AbsolutePoses_true(:,1:3,k)'*AbsolutePoses_true(:,4,k) + sigma_p*randn(3,1));
    AbsolutePoses(:,:,k) = [Ck_R_C1 Ck_t_C1];
end


sigma_fc = 30;
sigma_cc = 10;

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


usedPoses  = 1;

AbsolutePoses_EPnP = AbsolutePoses;
AbsolutePoses_EPnP_Gaussian = AbsolutePoses;
AbsolutePoses_EPnP_2 = AbsolutePoses;
AbsolutePoses_EPnP_Gaussian_2 = AbsolutePoses;
repro_error = zeros(NumOfPoses,5);
for k = 2:NumOfPoses
    InlierIdx_epnp = find(PoseGraphMatrix(:,k));
    x3d_h = FeatureBag_true_xyz(:, InlierIdx_epnp)';
    x2d_h = featureExtracted_pixel{k}';
    
    %EPnP
    [Rp,Tp,Xc,sol]=efficient_pnp(x3d_h,x2d_h,CameraParams.K);
    AbsolutePoses_EPnP(:,:,k) = [Rp Tp];
    repro_error(k,1)=reprojection_error_usingRT(x3d_h(:,1:3),x2d_h(:,1:2),Rp,Tp,CameraParams.K);
    
    %EPnP 2
    [Rp,Tp,Xc,sol]=efficient_pnp_2(x3d_h,x2d_h,CameraParams.K);
    AbsolutePoses_EPnP_2(:,:,k) = [Rp Tp];
    repro_error(k,4)=reprojection_error_usingRT(x3d_h(:,1:3),x2d_h(:,1:2),Rp,Tp,CameraParams.K);
    
    %EPnP w/ Gausssian 
    [Rp,Tp,Xc,sol]=efficient_pnp_gauss(x3d_h,x2d_h,CameraParams.K);
    AbsolutePoses_EPnP_Gaussian(:,:,k) = [Rp Tp];
    repro_error(k,2)=reprojection_error_usingRT(x3d_h(:,1:3),x2d_h(:,1:2),Rp,Tp,CameraParams.K);
    
    %EPnP w/ Gausssian 2
    [Rp,Tp,Xc,sol]=efficient_pnp_gauss_2(x3d_h,x2d_h,CameraParams.K);
    AbsolutePoses_EPnP_Gaussian_2(:,:,k) = [Rp Tp];
    repro_error(k,5)=reprojection_error_usingRT(x3d_h(:,1:3),x2d_h(:,1:2),Rp,Tp,CameraParams.K);
    
    % PnP
    [Cn_R_Cr, Cn_t_Cr, InlierIdx, ~, Cr_P, Cn_z,idx_3d] = simulation_P3P_RANSAC_InvDep_pixel(PoseGraphMatrix, AbsolutePoses,  ...
                                                         featureExtracted_pixel,                      ...
                                                         usedPoses, k,                      ...
                                                         FeatureBag_true_InvDepth, CameraParams,          ...
                                                         500, 0.7, 50); 
    [Cn_R_Cr, Cn_t_Cr, ~, normHist] = simulation_PnP_NL_InvDep_pixel(Cn_R_Cr, Cn_t_Cr,AbsolutePoses, FeatureBag_true_InvDepth(:, InlierIdx), featureExtracted_pixel{k}(1:2,InlierIdx),CameraParams);                  
    AbsolutePoses(:,:,k) = [Cn_R_Cr, Cn_t_Cr];
    repro_error(k,3)=reprojection_error_usingRT(FeatureBag_true_xyz(1:3, InlierIdx)',featureExtracted_pixel{k}(1:2,InlierIdx)',Cn_R_Cr,Cn_t_Cr,CameraParams.K);  

    usedPoses = [usedPoses k];
%     [FeatureBag_InvDep, AbsolutePoses_InvDep, CameraParams_InvDep,~, ~,~] = BA_Intrinsic_simulation(FeatureBag_true_InvDepth, AbsolutePoses, PoseGraphMatrix, featureExtracted_pixel, usedPoses, CameraParams);
end

% error_intrinsic = zeros(4,2);
% error_intrinsic(1,1) = abs(CameraParams_InvDep.fc(1) - K_true(1,1));
% error_intrinsic(2,1) = abs(CameraParams_InvDep.fc(2) - K_true(2,2));
% error_intrinsic(3,1) = abs(CameraParams_InvDep.cc(1) - K_true(1,3));
% error_intrinsic(4,1) = abs(CameraParams_InvDep.cc(2) - K_true(2,3));
% error_intrinsic(1,2) = abs(CameraParams.fc(1) - K_true(1,1));
% error_intrinsic(2,2) = abs(CameraParams.fc(2) - K_true(2,2));
% error_intrinsic(3,2) = abs(CameraParams.cc(1) - K_true(1,3));
% error_intrinsic(4,2) = abs(CameraParams.cc(2) - K_true(2,3));
% figure;
% plot(error_intrinsic(:,1),'r*');hold on
% plot(error_intrinsic(:,2),'b*');
% legend('output','input');
% title('error of intrinsic matrix');

% camera pose
figure;
error_r = zeros(NumOfPoses,4);
error_t = zeros(NumOfPoses,4);
for k = 1:NumOfPoses
    error_r(k,1) = acos((trace(AbsolutePoses_EPnP(:,1:3,k)*AbsolutePoses_true(:,1:3,k)')-1)/2);
    error_t(k,1) = sum((AbsolutePoses_EPnP(:,4,k)- AbsolutePoses_true(:,4,k)).^2)^0.5;
    error_r(k,2) = acos((trace(AbsolutePoses_EPnP_Gaussian(:,1:3,k)*AbsolutePoses_true(:,1:3,k)')-1)/2);
    error_t(k,2) = sum((AbsolutePoses_EPnP_Gaussian(:,4,k)- AbsolutePoses_true(:,4,k)).^2)^0.5;
    error_r(k,3) = acos((trace(AbsolutePoses(:,1:3,k)*AbsolutePoses_true(:,1:3,k)')-1)/2);
    error_t(k,3) = sum((AbsolutePoses(:,4,k)- AbsolutePoses_true(:,4,k)).^2)^0.5;
    error_r(k,4) = acos((trace(AbsolutePoses_EPnP_2(:,1:3,k)*AbsolutePoses_true(:,1:3,k)')-1)/2);
    error_t(k,4) = sum((AbsolutePoses_EPnP_2(:,4,k)- AbsolutePoses_true(:,4,k)).^2)^0.5;
    error_r(k,5) = acos((trace(AbsolutePoses_EPnP_Gaussian_2(:,1:3,k)*AbsolutePoses_true(:,1:3,k)')-1)/2);
    error_t(k,5) = sum((AbsolutePoses_EPnP_Gaussian_2(:,4,k)- AbsolutePoses_true(:,4,k)).^2)^0.5;
end
subplot(2,1,1);hold on
plot(error_r(:,1));hold on
plot(error_r(:,2));hold on
plot(error_r(:,3));hold on
plot(error_r(:,4));hold on
plot(error_r(:,5));
legend('EPnP','EPnP w/Gaussia','PnP','EPnP 2','EPnP w/Gaussia 2');
title('error of rotation');
subplot(2,1,2);hold on
plot(error_t(:,1));hold on
plot(error_t(:,2));hold on
plot(error_t(:,3));hold on
plot(error_t(:,4));hold on
plot(error_t(:,5));
legend('EPnP','EPnP w/Gaussian','PnP','EPnP 2','EPnP w/Gaussia 2');
title('error of position');

figure;
plot(repro_error(:,1));hold on
plot(repro_error(:,2));hold on
plot(repro_error(:,3));hold on
plot(repro_error(:,4));hold on
plot(repro_error(:,5));
legend('EPnP','EPnP w/Gaussia','PnP','EPnP 2','EPnP w/Gaussia 2');
title('reprojection error');


% reprojection
figure;
for i = 1:NumOfPoses
        subplot(NumOfPoses,1,i); hold on
        plot(featureExtracted_true_pixel{i}(1,:),featureExtracted_true_pixel{i}(2,:),'r*');hold on
%         axis([0 CameraParams.cc(1)*2 0 CameraParams.cc(2)*2]);hold on
        triangulate_feat = find(FeatureBag_true_xyz(4,:) & PoseGraphMatrix(:,i)');
        repro = zeros(3,length(triangulate_feat));
        for k = triangulate_feat
            % PnP
%             m_theta_phi = [cos(FeaturesBag_InvDep(2,k))*cos(FeaturesBag_InvDep(1,k));cos(FeaturesBag_InvDep(2,k))*sin(FeaturesBag_InvDep(1,k));sin(FeaturesBag_InvDep(2,k))];
            Ck_P_W = AbsolutePoses(:,:,i)*[InversePose(AbsolutePoses(:,:,FeatureBag_true_xyz(4,k)));0 0 0 1]*[FeatureBag_true_xyz(1:3,k);1];
            Ck_P_W = Ck_P_W/Ck_P_W(3);
            repro = CameraParams.K*Ck_P_W;
            plot(repro(1),repro(2),'bo');hold on
            % EPnP
            Ck_P_W = AbsolutePoses_EPnP(:,:,i)*[InversePose(AbsolutePoses_EPnP(:,:,FeatureBag_true_xyz(4,k)));0 0 0 1]*[FeatureBag_true_xyz(1:3,k);1];
            Ck_P_W = Ck_P_W/Ck_P_W(3);
            repro = CameraParams.K*Ck_P_W;
            plot(repro(1),repro(2),'mo');hold on
            % EPnP w/ Gaussian
            Ck_P_W = AbsolutePoses_EPnP_Gaussian(:,:,i)*[InversePose(AbsolutePoses_EPnP_Gaussian(:,:,FeatureBag_true_xyz(4,k)));0 0 0 1]*[FeatureBag_true_xyz(1:3,k);1];
            Ck_P_W = Ck_P_W/Ck_P_W(3);
            repro = CameraParams.K*Ck_P_W;
            plot(repro(1),repro(2),'go');hold on
        end
        legend('measurement','PnP','EPnP','EPnP w/ Gaussian');
        title(['reprojection of Image ',num2str(i)]);
end
