function repro_error = ComputeReprojectionError(PoseGraphMatrix,AbsolutePoses,FeaturesBag,featureExtracted,feat_view_3d_idx,CameraParams)

global CameraModel
% feat_view_3d_idx = 95;
cam_idx = find(PoseGraphMatrix(feat_view_3d_idx,:));
repro_error = zeros(length(cam_idx),1);

% figure;
for newPose = 1:length(cam_idx)
    repro_pose = cam_idx(newPose);
    
    feat_view_2d_idx = PoseGraphMatrix(feat_view_3d_idx,repro_pose);
    Ck_T_C1 = AbsolutePoses(:,:,repro_pose);
    Cr_T_C1 = AbsolutePoses(:,:,FeaturesBag(4,feat_view_3d_idx));
    Ck_T_Cr = Ck_T_C1*[InversePose(Cr_T_C1);zeros(1,3) 1];
    feat_2d = featureExtracted{repro_pose}(:,feat_view_2d_idx);
    
    
    F_xyz = 1/FeaturesBag(3,feat_view_3d_idx)*[cos(FeaturesBag(2,feat_view_3d_idx))*cos(FeaturesBag(1,feat_view_3d_idx));...
             cos(FeaturesBag(2,feat_view_3d_idx))*sin(FeaturesBag(1,feat_view_3d_idx));sin(FeaturesBag(2,feat_view_3d_idx))];
    feat_repro = Ck_T_Cr*[F_xyz;1];
    feat_repro(1:2,:) = feat_repro(1:2,:)./feat_repro(3,:);
    feat_repro = StaticDistort(CameraParams.fc,CameraParams.cc,CameraParams.kc,feat_repro(1:2,:),CameraModel);
    repro_error(newPose,:) = norm(feat_repro(1:2,:) - feat_2d(1:2,:));
    
%     subplot(length(cam_idx),1,newPose);
%     plot(feat_repro(1,:),feat_repro(2,:),'bo');hold on
%     plot(feat_2d(1,:),feat_2d(2,:),'r*');
%     legend('estimate','measurement');
%     title(['camera ',num2str(repro_pose)]);
end