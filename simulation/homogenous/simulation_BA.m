clear all; close all;clc
codepath = '/home/marsfei/Fei/PTAM';

addpath(codepath);
addpath([codepath,'/Drawings/']);
addpath([codepath,'/robotics3D']);
addpath([codepath,'/JacobianModule/']);
addpath([codepath,'/BundleAdjustment']);
addpath([codepath,'/RobustCostFunction/']);
addpath([codepath,'/InverseDepth/']);
addpath([codepath,'/TwoViewReconstruction/']);

NumOfFeatures = 100;
NumOfPoses = 8;


% PoseGraphMatrix = ones(NumOfFeatures, NumOfPoses);
PoseGraphMatrix = randi([0 1], NumOfFeatures, NumOfPoses);
mintracklength = 4;
addMoreFeature = find(sum(PoseGraphMatrix,2)<=mintracklength);
% Make sure each feature view by at least 4 images
for k = 1:length(addMoreFeature)
    cam_not_view = find(PoseGraphMatrix(addMoreFeature(k),:)==0);
    count = length(find(PoseGraphMatrix(addMoreFeature(k),:)));
    for i = 1:mintracklength-count
        PoseGraphMatrix(addMoreFeature(k),cam_not_view(i)) = 1;
    end
end



%% 3D-WORLD creation


% To ensure f.o.v = 60 deg.
fov = 60*pi/180;
Rfeat = 1; Rcam = 3*Rfeat/sin(fov/2);


% Random pose generator
% W_T_C1 = [];
AbsolutePoses_true = zeros(3,4,NumOfPoses);
for k = 1:NumOfPoses
    tt = pi/2*(k-1)/NumOfPoses;
    WRC = [-sin(tt) 0 -cos(tt); cos(tt) 0 -sin(tt); 0 -1 0];
    WpC = Rcam*[cos(tt);sin(tt);0];
    
%     AbsolutePoses_true(:,:,k) = [WRC WpC];
    
    
    if k == 1
        W_T_C1 = [WRC WpC];
        AbsolutePoses_true(:,:,k) = [eye(3) zeros(3,1)];
    else
        C_T_W = [WRC' -WRC'*WpC];
        AbsolutePoses_true(:,:,k) = C_T_W*[W_T_C1;zeros(1,3) 1];
    end
end



% Random feature generator
FeatureBag_true = zeros(4,NumOfFeatures);
FeatureBag_true_xyz = zeros(4,NumOfFeatures);
%figure(2); hold on;
for k = 1:NumOfFeatures
    tt = 2*pi*rand();
    phi = pi*rand() - pi/2;
    
    First_view_cam_idx = min(find(PoseGraphMatrix(k,:))); % first camera view this feature
%     if isempty(First_view_cam_idx)
%         First_view_cam_idx = 0;
%     end
%     FeatureBag_true(:,k) = [tt;phi;Rfeat;First_view_cam_idx];
    FeatureBag_true_xyz(:,k) = [Rfeat*[cos(phi)*cos(tt);cos(phi)*sin(tt);sin(phi)];First_view_cam_idx]; % <--- w.r.t the world
    FeatureBag_true_xyz(1:3,k) = AbsolutePoses_true(:,:,First_view_cam_idx)*[InversePose(W_T_C1);zeros(1,3) 1]*[FeatureBag_true_xyz(1:3,k);1];
    
    % sphere coordinate
    FeatureBag_true(3,k) = 1/norm(FeatureBag_true_xyz(1:3,k));
    FeatureBag_true(1,k) = atan2(FeatureBag_true_xyz(2,k),FeatureBag_true_xyz(1,k));
    FeatureBag_true(2,k) = atan2(FeatureBag_true_xyz(3,k),(FeatureBag_true_xyz(1,k)^2+FeatureBag_true_xyz(2,k)^2)^0.5);
    FeatureBag_true(4,k) = First_view_cam_idx;   
end



VisualizeMultiPoses(AbsolutePoses_true, FeatureBag_true_xyz, 1:(NumOfPoses-1), NumOfPoses);


% Create feature projections as measurement
featureExtracted_true = cell(NumOfPoses, 1);
for k = 1:NumOfPoses    
    numFeatPosek = find(PoseGraphMatrix(:,k));
    featureExtracted_true{k} = zeros(2, length(numFeatPosek));
    for l = 1:length(numFeatPosek)
        % assign correct feature index for pose graph matrix
        PoseGraphMatrix(numFeatPosek(l),k) = l; 
        
        % Pose Transformation
        Ck_T_W = AbsolutePoses_true(:,:,k);
        Cr_T_W = AbsolutePoses_true(:,:,FeatureBag_true_xyz(4,numFeatPosek(l)));
        Ck_T_Cr = Ck_T_W*[InversePose(Cr_T_W);zeros(1,3) 1];
        
        % Image Projection
        Ck_p_f = Ck_T_Cr*[FeatureBag_true_xyz(1:3,numFeatPosek(l));1];
        featureExtracted_true{k}(:,l) = Ck_p_f(1:2)./Ck_p_f(3);
    end
    featureExtracted_true{k} = [featureExtracted_true{k};ones(1,size(featureExtracted_true{k},2))];
end

%% add noise 
sigma_z = 0; % noise for measurement
sigma_theta = 0.02; % 5e-2, noise for 3d feature 
sigma_phi = 0.02; % 5e-2, noise for 3d feature
sigma_rho = 0.01; % 5e-2, noise for 3d feature
sigma_q = 0.01; % noise for current camera orientation
sigma_p = 0.1; % 0.2, noise for current camera position

featureExtracted = cell(NumOfPoses, 1);
f_tilde = cell(NumOfPoses, 1);
for k = 1:NumOfPoses
   f_tilde{k} = sigma_z*randn(size(featureExtracted_true{k}));
   featureExtracted{k} = featureExtracted_true{k} + f_tilde{k};
end

FeatureBag = zeros(4,NumOfFeatures);
FeatureBag_xyz = zeros(4,NumOfFeatures);
F_tilde = zeros(3,NumOfFeatures);
for k = 1:NumOfFeatures
    F_tilde(:,k) = randn(3,1).*[sigma_theta;sigma_phi;sigma_rho];
    FeatureBag(1:3,k) = FeatureBag_true(1:3,k) + F_tilde(:,k);   
    FeatureBag(4,k) = FeatureBag_true(4,k);
    FeatureBag_xyz(1:3,k) = 1/FeatureBag(3,k)*[cos(FeatureBag(2,k))*cos(FeatureBag(1,k));cos(FeatureBag(2,k))*sin(FeatureBag(1,k));sin(FeatureBag(2,k))];
    FeatureBag_xyz(4,k) = FeatureBag(4,k);
end

seen_feat = cell(NumOfPoses,1);
for k = 1:NumOfPoses
    seen_feat{k} = find(PoseGraphMatrix(:,k));
end

AbsolutePoses = zeros(3,4,NumOfPoses);
q_tilde = zeros(3,NumOfPoses);
p_tilde = zeros(3,NumOfPoses);


for k = 1
    AbsolutePoses(:,:,k) = AbsolutePoses_true(:,:,k);
end

for k = 2:NumOfPoses
    q_tilde(:,k) = randn(3,1)*sigma_q;
    dq = [0.5*q_tilde(:,k);1];
    dq = dq / norm(dq);   
    p_tilde(:,k) = sigma_p*randn(3,1);
    Ck_R_W = quat2rot(quat_mul(dq, rot2quat(AbsolutePoses_true(:,1:3,k))));
    Ck_p_W = AbsolutePoses_true(:,4,k) - AbsolutePoses_true(:,1:3,k)*p_tilde(:,k);
    AbsolutePoses(:,:,k) = [Ck_R_W Ck_p_W];
end

pertube.q_tilde = q_tilde;
pertube.p_tilde = p_tilde;
pertube.F_tilde = F_tilde;


CameraParams.Kinv = eye(3);
CameraParams.K = eye(3);
usedPoses = 1:NumOfPoses;
[FeatureBag_xyz, AbsolutePoses_xyz, ~, ~,~] = BA_Cartesian_2viewBLS_simulation(FeatureBag_xyz, AbsolutePoses, PoseGraphMatrix, featureExtracted, usedPoses, CameraParams);
[FeatureBag_xyz, AbsolutePoses_xyz, ~, ~,~] = BA_Cartesian_2viewBLS_Robust(FeatureBag_xyz, AbsolutePoses, PoseGraphMatrix, featureExtracted, usedPoses, CameraParams, 1e-1);

% [FeatureBag_InvDep, AbsolutePoses_InvDep, ~, ~,~,error] = BA_InverseDepth_simulation(FeatureBag, AbsolutePoses, PoseGraphMatrix, featureExtracted, usedPoses, CameraParams);

error_r = zeros(NumOfPoses,3);
error_t = zeros(NumOfPoses,3);
for k = 1:NumOfPoses
%     AbsolutePoses_xyz(:,1:3,k) = AbsolutePoses_xyz(:,1:3,k)';
%     AbsolutePoses_xyz(:,4,k) = -AbsolutePoses_xyz(:,1:3,k)*AbsolutePoses_xyz(:,4,k);
%     AbsolutePoses_InvDep(:,1:3,k) = AbsolutePoses_InvDep(:,1:3,k)';
%     AbsolutePoses_InvDep(:,4,k) = -AbsolutePoses_InvDep(:,1:3,k)*AbsolutePoses_InvDep(:,4,k);
    
    error_r(k,1) = acos((trace(AbsolutePoses_xyz(:,1:3,k)*AbsolutePoses_true(:,1:3,k)')-1)/2);
    error_t(k,1) = sum((AbsolutePoses_xyz(:,4,k)- AbsolutePoses_true(:,4,k)).^2)^0.5;
    error_r(k,2) = acos((trace(AbsolutePoses_InvDep(:,1:3,k)*AbsolutePoses_true(:,1:3,k)')-1)/2);
    error_t(k,2) = sum((AbsolutePoses_InvDep(:,4,k)- AbsolutePoses_true(:,4,k)).^2)^0.5;
    
    error_r(k,3) = acos((trace(AbsolutePoses(:,1:3,k)*AbsolutePoses_true(:,1:3,k)')-1)/2);
    error_t(k,3) = sum((AbsolutePoses(:,4,k)- AbsolutePoses_true(:,4,k)).^2)^0.5;
end
% 
figure;
subplot(2,1,1);
plot(error_r(:,1));hold on
plot(error_r(:,2));hold on
plot(error_r(:,3));
legend('Cartisan','Inverse Depth','Input');
title('rotation');
subplot(2,1,2);
plot(error_t(:,1));hold on
plot(error_t(:,2));hold on
plot(error_t(:,3));
legend('Cartisan','Inverse Depth','Input');
title('position');
