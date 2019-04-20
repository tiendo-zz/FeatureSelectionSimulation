function [FeatureBag_true_InvDepth, FeatureBag_true_xyz, featureExtracted_true, AbsolutePoses_true, PoseGraphMatrix] = GenerateGroundTruth_ball(NumOfFeatures, NumOfPoses)

PoseGraphMatrix = ones(NumOfFeatures, NumOfPoses);
% PoseGraphMatrix = randi([0 1], NumOfFeatures, NumOfPoses);
% mintracklength = 6;
% addMoreFeature = find(sum(PoseGraphMatrix,2)<=mintracklength);
% % Make sure each feature view by at least 4 images
% for k = 1:length(addMoreFeature)
%     cam_not_view = find(PoseGraphMatrix(addMoreFeature(k),:)==0);
%     count = length(find(PoseGraphMatrix(addMoreFeature(k),:)));
%     for i = 1:mintracklength-count
%         PoseGraphMatrix(addMoreFeature(k),cam_not_view(i)) = 1;
%     end
% end



%% 3D-WORLD creation


% To ensure f.o.v = 60 deg.
fov = 60*pi/180;
Rfeat = 1; Rcam = 3*Rfeat/sin(fov/2);


% Random pose generator
% W_T_C1 = [];
AbsolutePoses_true = zeros(3,4,NumOfPoses);
for k = 1:NumOfPoses
    tt = pi*(k-1)/NumOfPoses;
    R_z = [cos(tt) -sin(tt) 0;sin(tt) cos(tt) 0;0 0 1];
    WRC = [-sin(tt) 0 -cos(tt); cos(tt) 0 -sin(tt); 0 -1 0]*R_z;
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
FeatureBag_true_InvDepth = zeros(4,NumOfFeatures);
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
    FeatureBag_true_InvDepth(3,k) = 1/norm(FeatureBag_true_xyz(1:3,k));
    FeatureBag_true_InvDepth(1,k) = atan2(FeatureBag_true_xyz(2,k),FeatureBag_true_xyz(1,k));
    FeatureBag_true_InvDepth(2,k) = atan2(FeatureBag_true_xyz(3,k),(FeatureBag_true_xyz(1,k)^2+FeatureBag_true_xyz(2,k)^2)^0.5);
    FeatureBag_true_InvDepth(4,k) = First_view_cam_idx;   
end


% VisualizeMultiPoses(AbsolutePoses_true, FeatureBag_true_xyz, 1:(NumOfPoses-1), NumOfPoses);


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