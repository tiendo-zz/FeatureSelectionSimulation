function [FeatureBag_true_InvDepth, FeatureBag_true_xyz, featureExtracted_true, AbsolutePoses_true, PoseGraphMatrix] = GenerateGroundTruth_Circle(NumOfFeatures, NumOfPoses)

x_range = 10;
y_range = 10;

% generate camera pose

AbsolutePoses_true = zeros(3,4,NumOfPoses);
alpha = 1;
for k = 1:NumOfPoses
    WRC = [-1 0 0;0 1 0;0 0 -1];
    theta = k/NumOfPoses*2*pi;
    WpC = [x_range/2+x_range/4*cos(theta);y_range/2+y_range/4*sin(theta);10-theta*alpha];
    
%     AbsolutePoses_true(:,:,k) = [WRC' -WRC'*WpC];
    if k == 1
        W_T_C1 = [WRC WpC];
        AbsolutePoses_true(:,:,k) = [eye(3) zeros(3,1)];
    else
        C_T_W = [WRC' -WRC'*WpC];
        AbsolutePoses_true(:,:,k) = C_T_W*[W_T_C1;zeros(1,3) 1];
    end
end

% generate 3d features
FeatureBag_true_W = zeros(4,NumOfFeatures);
FeatureBag_true_W(1:3, :) = [x_range*rand(1,NumOfFeatures);...
                               y_range*rand(1,NumOfFeatures);...
                               zeros(1,NumOfFeatures)];
FeatureBag_true_W(4, :) = ones(1,NumOfFeatures);
FeatureBag_true_xyz = zeros(4,NumOfFeatures);
FeatureBag_true_xyz(1:3,:) = InversePose(W_T_C1)*FeatureBag_true_W;
                               


PoseGraphMatrix = zeros(NumOfFeatures,NumOfPoses);


% create 2d-3d correspondence
fov = 90*pi/180;
FeatureBag_true_InvDepth = zeros(4,NumOfFeatures);
for i = 1:NumOfFeatures
    flag = 0;
    for k = 1:NumOfPoses
        Ck_p_f = AbsolutePoses_true(:,:,k)*[FeatureBag_true_xyz(1:3,i);1];
        phi = atan2(Ck_p_f(3),(Ck_p_f(1)^2+Ck_p_f(2)^2)^0.5);
        if (pi/2 - phi) < fov/2
            if flag == 0
               flag = 1;
               First_view_cam_idx = k;
               FeatureBag_true_xyz(:,i) = [Ck_p_f(1:3);First_view_cam_idx];
               FeatureBag_true_InvDepth(3,i) = 1/norm(Ck_p_f(1:3));
               FeatureBag_true_InvDepth(1,i) = atan2(Ck_p_f(2),Ck_p_f(1));
               FeatureBag_true_InvDepth(2,i) = atan2(Ck_p_f(3),(Ck_p_f(1)^2+Ck_p_f(2)^2)^0.5);
               FeatureBag_true_InvDepth(4,i) = First_view_cam_idx;               
            end  
            PoseGraphMatrix(i,k) = 1;
        end        
    end
end
        

% VisualizeMultiPoses(AbsolutePoses_true, FeatureBag_true_xyz, 1:(NumOfPoses-1), NumOfPoses);



% create 2d features
featureExtracted_true = cell(NumOfPoses, 1);
for k = 1:NumOfPoses    
    numFeatPosek = find(PoseGraphMatrix(:,k));
    featureExtracted_true{k} = zeros(2, length(numFeatPosek));
    for l = 1:length(numFeatPosek)
        % assign correct feature index for pose graph matrix
        PoseGraphMatrix(numFeatPosek(l),k) = l; 
        
        % Pose Transformation
        Ck_T_C1 = AbsolutePoses_true(:,:,k);
        Cr_T_C1 = AbsolutePoses_true(:,:,FeatureBag_true_xyz(4,numFeatPosek(l)));
        Ck_T_Cr = Ck_T_C1*[InversePose(Cr_T_C1);zeros(1,3) 1];
        
        % Image Projection
        Ck_p_f = Ck_T_Cr*[FeatureBag_true_xyz(1:3,numFeatPosek(l));1];
        featureExtracted_true{k}(:,l) = Ck_p_f(1:2)./Ck_p_f(3);
    end
    featureExtracted_true{k} = [featureExtracted_true{k};ones(1,size(featureExtracted_true{k},2))];
end