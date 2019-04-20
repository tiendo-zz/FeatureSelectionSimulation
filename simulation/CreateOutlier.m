function [featureExtracted_outlier,PoseGraphMatrix_outlier] = CreateOutlier(PoseGraphMatrix,featureExtracted,outlier_num,noise_level)

observe_measurement = (PoseGraphMatrix~=0);
featureExtracted_outlier = cell(size(PoseGraphMatrix,2),1);
PoseGraphMatrix_outlier = zeros(size(PoseGraphMatrix));

outlier_num_per_img = round(outlier_num/size(PoseGraphMatrix,2));
for pose_id = 1:size(PoseGraphMatrix,2)

    observe_pose_k = find(observe_measurement(:,pose_id));
    
    if  length(observe_pose_k) < outlier_num_per_img
        featureExtracted_outlier{pose_id} = featureExtracted{pose_id};
        continue;
    end
    
    outlier_landmark_id = observe_pose_k(randperm(length(observe_pose_k),outlier_num_per_img));
    featureExtracted_outlier{pose_id} = featureExtracted{pose_id};
    
    for outlier = outlier_landmark_id'
        observation_id = PoseGraphMatrix(outlier,pose_id);
        PoseGraphMatrix_outlier(outlier,pose_id) = observation_id;
        featureExtracted_outlier{pose_id}(:,observation_id) = [noise_level + featureExtracted{pose_id}(1:2,observation_id);1];
    end
end