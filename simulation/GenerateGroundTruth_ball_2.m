function [FeatureBag_true_InvDepth, FeatureBag_true_xyz, featureExtracted_true, AbsolutePoses_true, PoseGraphMatrix] = GenerateGroundTruth_ball_2(angle_range,NumOfFeatures, NumOfPoses, K)

PoseGraphMatrix = zeros(NumOfFeatures, NumOfPoses);
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
% quickly compute the max fov %% TODO: use closed form formula:
x_dim = 2 * K(1,3); y_dim = 2 * K(2,3);
sample_homo = K \ [0;0;1];
fov = abs(atan2( norm(sample_homo(1:2)), 1) );
Rfeat = 30; Rcam = Rfeat/sin(fov);


% Random pose generator
% W_T_C1 = [];
AbsolutePoses_true = zeros(3,4,NumOfPoses);
for k = 1:NumOfPoses
    tt = angle_range*(k-1)/NumOfPoses;
    phi = angle_range*(k-1)/NumOfPoses;
    % WRC = [-sin(tt) 0 -cos(tt); cos(tt) 0 -sin(tt); 0 -1 0];
    WpC = Rcam*[cos(tt)*cos(phi);sin(tt)*cos(phi);sin(phi)];
    k_axis = skewsymm(WpC)* [0;0;-1]; k_axis = k_axis/norm(k_axis);
    CRW = aa2rot(k_axis, acos(-WpC(3)/norm(WpC)));
    WRC = CRW';
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
k = 1;
while( k <= NumOfFeatures )
    PoseGraphMatrix(k, :) = zeros(1, NumOfPoses);
    tt = 2*pi*rand();
    phi = pi*rand() - pi/2;
    W_p_f = (0.5 + 0.5*rand())*Rfeat*[cos(phi)*cos(tt);cos(phi)*sin(tt);sin(phi)]; % <--- w.r.t the world    
    flag = 0;
    
    observation_count = 0;
    for j = 1:NumOfPoses
        Cj_p_f = AbsolutePoses_true(:,:,j)*[InversePose(W_T_C1);zeros(1,3) 1]*[W_p_f;1];        
        pixel_value = K * (Cj_p_f / Cj_p_f(3));
        if pixel_value(1) <= x_dim && pixel_value(1) >= 0 && pixel_value(2) <= y_dim && pixel_value(2) >= 0
            if flag == 0
                flag = 1;
                First_view_cam_idx = j;
                FeatureBag_true_xyz(:,k) = [Cj_p_f;First_view_cam_idx];
                FeatureBag_true_InvDepth(3,k) = 1/norm(Cj_p_f(1:3));
                FeatureBag_true_InvDepth(1,k) = atan2(Cj_p_f(2),Cj_p_f(1));
                FeatureBag_true_InvDepth(2,k) = atan2(Cj_p_f(3),(Cj_p_f(1)^2+Cj_p_f(2)^2)^0.5);
                FeatureBag_true_InvDepth(4,k) = First_view_cam_idx; 
            end
            PoseGraphMatrix(k, j) = 1;    
            observation_count = observation_count + 1;
        end
    end
    if observation_count == NumOfPoses
        k = k + 1; % Next feature generating
    end
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