function   [FeaturesBag,PoseGraphMatrix] = MultiPoseTriangulation_InvDep(FeaturesBag,           ...
                                                 AbsolutePoses,         ...
                                                 PoseGraphMatrix,       ...
                                                 newPose, usedPoses,    ...
                                                 featureExtracted, ...                                         
                                                 CameraParams,threshold)
                                             
global CameraModel                                             
                                             
newposeIdx = find(PoseGraphMatrix(:, newPose));
usedposeIdx = find(sum(PoseGraphMatrix(:, usedPoses),2));
triangulateIdx = intersect(newposeIdx, usedposeIdx);
% triangulateIdx = setdiff(triangulateIdx,find(FeaturesBag(4,:)~=0));


% MultiPoses linear triangulation
allPoses = [usedPoses, newPose];
for k = triangulateIdx'
    % 1. Find how many poses see this feature
    % (need at least three poses viewing this feature)
    posesViewkfeatIdx = find(PoseGraphMatrix(k,allPoses));    

    [temp_feat,temp_feat_xyz] = LinearMultiPoseTriangulation_InvDep(k,AbsolutePoses,         ...
                                                 allPoses, ...
                                                 PoseGraphMatrix,       ...
                                                 featureExtracted, ...                                         
                                                 CameraParams,CameraModel);
    if isempty(temp_feat_xyz)
        continue;
    end

    if( ChieralityCheckMultiPoses(temp_feat_xyz, AbsolutePoses, allPoses(posesViewkfeatIdx)) )
       triangulate_decision = ReprojectionMultiPose(k,PoseGraphMatrix,...
                                                    allPoses,...
                                                    AbsolutePoses,...
                                                    featureExtracted,...
                                                    temp_feat_xyz,...
                                                    CameraParams,...]
                                                    CameraModel,threshold);
        if triangulate_decision == 1
            FeaturesBag(:,k) = temp_feat;
        else
            % if current feature is triangulated before
            if(FeaturesBag(4,k) ~=0) 
               % erase measurement of new pose
               PoseGraphMatrix(k,newPose) = 0;
               % erase the measurement of previous pose added because of new
               % pose
               for pose_view_feat = find(PoseGraphMatrix(k,:))
                   if pose_view_feat < FeaturesBag(4,k)
                      PoseGraphMatrix(k,pose_view_feat) = 0;
                   end
               end
               FeaturesBag(:,k) = 0;
               [temp_feat,temp_feat_xyz] = LinearMultiPoseTriangulation_InvDep(k,AbsolutePoses,         ...
                                                 allPoses, ...
                                                 PoseGraphMatrix,       ...
                                                 featureExtracted, ...                                         
                                                 CameraParams,CameraModel);
               if ~isempty(temp_feat_xyz)
                   if( ChieralityCheckMultiPoses(temp_feat_xyz, AbsolutePoses, allPoses(posesViewkfeatIdx)) )
                       triangulate_decision = ReprojectionMultiPose(k,PoseGraphMatrix,...
                                                                    allPoses,...
                                                                    AbsolutePoses,...
                                                                    featureExtracted,...
                                                                    temp_feat_xyz,...
                                                                    CameraParams,...]
                                                                    CameraModel,threshold);
                        if triangulate_decision == 1
                           FeaturesBag(:,k) = temp_feat;
                        end
                   end
               end           
            end
        end
    end
end


