function     [Rbest, tbest, bestInlier, bestGeneratorSet, featureUsed, measurementUsed,idx_3d] = P3P_RANSAC_InvDep_pixel(PoseGraphMatrix, ...
                                                                           AbsolutePoses,    ...
                                                                           featureExtracted, ...
                                                                           usedPoses, newPose,         ...
                                                                           FeaturesBag, CameraParams, ...
                                                                           MaxIter, inlier_ratio,threshold,P3P_inlier_thresh)  
% 2D-3D correspondences
% extract the indices of features seen in refPose in posegraph matrix that
% is already triangulated
% Rbest, tbest : w.r.t camera

global CameraModel

P = zeros(3,1000); 
P_InvDep = zeros(4,1000); 
w = zeros(1000,2);
idx_3d = zeros(1000,1);
count = 0;
for posei = usedPoses
    available_features_idx = find(FeaturesBag(4,:) == posei);
    if ~isempty(available_features_idx)
      [w_full, common_feature_indices] = GetMatchFromPoseGraph(PoseGraphMatrix, available_features_idx, posei, newPose, ...
                               featureExtracted{posei}, featureExtracted{newPose});
      
      % transform to Cartisian coordinate
      F_xyz = zeros(3,length(common_feature_indices));
      for k = 1:length(common_feature_indices)
          F_xyz(:,k) =  1/FeaturesBag(3,common_feature_indices(k))*[cos(FeaturesBag(2,common_feature_indices(k)))*cos(FeaturesBag(1,common_feature_indices(k)));...
              cos(FeaturesBag(2,common_feature_indices(k)))*sin(FeaturesBag(1,common_feature_indices(k)));sin(FeaturesBag(2,common_feature_indices(k)))];
      end
                           
      P(:,count+1:count+length(common_feature_indices)) = AbsolutePoses(:,1:3,posei)'*(F_xyz - ...
                    repmat(AbsolutePoses(:,4,posei), [1 length(common_feature_indices)]));
       P_InvDep(:,count+1:count+length(common_feature_indices)) = FeaturesBag(:,common_feature_indices);
       
       w(count+1:count+length(common_feature_indices),:) = w_full(:,4:5);
       idx_3d(count+1:count+length(common_feature_indices),:) = common_feature_indices;
       count = count + length(common_feature_indices);
    end
end
P(:,count+1:end) = [];
P_InvDep(:,count+1:end) = [];
w(count+1:end,:) = [];
idx_3d(count+1:end,:) = [];


% fprintf('P3P corresponding : %f\n', length(idx_3d)); 

featureUsed = P_InvDep;
measurementUsed = w;


% RANSAC parameter

MaxInlier = inlier_ratio*size(w,1);
iter = 1;
bestInlier = [];
if size(w,1) >= P3P_inlier_thresh
dist_thresh = 20;
    
while iter < MaxIter
    
    % Check uniqueness
    % Picking 4 distinct points
    picked_points_idx = randperm(size(w,1),4);
    flag = 1;
    for first_point = picked_points_idx
        for second_point = setdiff(picked_points_idx,first_point)
            dist = norm(w(first_point,:) - w(second_point,:));
            if dist < dist_thresh
                flag = 0;
                break;
            end
        end
        if flag == 0 
            break;
        end
    end
    if flag == 0
        continue;
    end
    
    % Compute pose  
    [C,p] = P3P_Solver(P(1:3,picked_points_idx(1:3)), StaticUndistort(CameraParams.fc,CameraParams.cc,CameraParams.kc,w(picked_points_idx(1:3), :)',CameraModel));
    if isempty(C) || isempty(p)
        continue;
    end
    [R,t,~]=TruePose_P3P_pixel(C,p,P(1:3,picked_points_idx(4)),w(picked_points_idx(4), :)',CameraParams,CameraModel);
    
    % Compute inlier based on reprojection
    [inlierIdx,error_mat] = ComputeInlier3D2DReprojection_pixel(P, w, CameraParams, R, t, threshold,CameraModel);
    
    if length(inlierIdx) > MaxInlier
        Rbest = R;
        tbest = t;
        bestInlier = inlierIdx;
        bestGeneratorSet = picked_points_idx;
        break;
    else
        if length(inlierIdx) > length(bestInlier)
            Rbest = R;
            tbest = t;
            bestInlier = inlierIdx;
            bestGeneratorSet = picked_points_idx;
        end
    end 
    
    iter = iter + 1;
end
else
    Rbest = []; tbest = [];
    bestGeneratorSet = [];
    bestInlier = [];
end
    