function [PoseGraphMatrix, pointerToImages,feature_view_b4_delete] = ConstructPoseGraphMatrix_p3p(InlierIdx, idx_3d, PoseGraphMatrix, pointerToImages, img_idx, FeatureBag)

% update pose graph matrix after p3p

% 1. extract outlier 3d idx
outlier_3d_idx = idx_3d(setdiff(1:length(idx_3d),InlierIdx));

% 2. if there are more than 4 imgs viewed this feature, clear current
% images' feature only; otherwise, delete this feature

count_1 = [];
count_2 = [];
if ~isempty(outlier_3d_idx)
   for m = 1:length(outlier_3d_idx)
       if sum(PoseGraphMatrix(outlier_3d_idx(m),:)~=0) >= 3
          % erase the measurement of previous pose added because of new
          % pose
          for pose_view_feat = find(PoseGraphMatrix(outlier_3d_idx(m),:))
              if pose_view_feat < FeatureBag(4,outlier_3d_idx(m))
                 PoseGraphMatrix(outlier_3d_idx(m),pose_view_feat) = 0;
              end
          end
          % the index in new image needed to be delete
          count_1 = [count_1;outlier_3d_idx(m)];
       else
          % the point in all image needed to be delete
          count_2 = [count_2;outlier_3d_idx(m)];
       end
   end
   if ~isempty(count_1)
      PoseGraphMatrix(count_1,img_idx) = 0;
   end
   if ~isempty(count_2)
       % Since we need to eliminate some entire rows, we need to
       % delete the numbers in pointerToImages.
       PoseGraphMatrix(count_2,:) = 0;
   end                 
end


% Update pointers
mess_row_pre = find(sum(PoseGraphMatrix(1:pointerToImages(img_idx),:),2)==0);
mess_row = find(sum(PoseGraphMatrix,2)==0);

feature_view_b4_delete = mess_row;
pointerToImages(img_idx) =  pointerToImages(img_idx) - length(mess_row_pre);
pointerToImages(img_idx+1) =  pointerToImages(img_idx+1) - length(mess_row);
% clear pose graph matrix
PoseGraphMatrix(mess_row,:) = [];
PoseGraphMatrix(pointerToImages(img_idx+1)+1:end, :) = []; 

