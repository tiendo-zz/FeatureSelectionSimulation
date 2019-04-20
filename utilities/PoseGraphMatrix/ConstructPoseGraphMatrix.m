function [PoseGraphMatrix, pointerToImages,feature_view_b4_delete] = ConstructPoseGraphMatrix(featureExtracted, PairWiseMatches, PoseGraphMatrix_old, pointerToImages_old, img_idx,match_image_id,flag)

% flag : 1 delete mess rows


PoseGraphMatrix = zeros(size(PoseGraphMatrix_old,1)+size(featureExtracted{img_idx},2), img_idx, 'uint16');
pointerToImages = zeros(img_idx+1,1);
% copy the old one to the new one
PoseGraphMatrix(1:size(PoseGraphMatrix_old,1),1:size(PoseGraphMatrix_old,2)) = PoseGraphMatrix_old;
pointerToImages(1:size(pointerToImages_old,1),1:size(pointerToImages_old,2)) = pointerToImages_old;
% pointer_reduce = 0;

% 1. Check if current image is the first image
if img_idx == 1
    new_features_Imagei = 1:size(featureExtracted{img_idx},2);
    PoseGraphMatrix(1+pointerToImages(img_idx):size(featureExtracted{img_idx},2)+pointerToImages(img_idx),img_idx) = ...
            (1:size(featureExtracted{img_idx},2))';
else
    for i = 1:length(match_image_id)
        matches_i = PairWiseMatches{SUTrilInd(img_idx,match_image_id(i),img_idx)};
        if isempty(matches_i)
            continue;
        end
        % 2. find other images who has the same features with current image,
        %    add the corresponding 2d feature id
       [idx_ij_of_new_features_a, idx_ij_of_new_features_b] = ismember(matches_i(:,3), PoseGraphMatrix(:,match_image_id(i)));
       % check if there are same preview features matching with different new
       % features, set this row to be zero.
       idx_ij_of_new_features_b = idx_ij_of_new_features_b(idx_ij_of_new_features_b~=0);
       idx_ij_of_new_features_a = find(idx_ij_of_new_features_a~=0);
       
       missMatch = find((PoseGraphMatrix(idx_ij_of_new_features_b, img_idx)~=0)&(PoseGraphMatrix(idx_ij_of_new_features_b, img_idx)~=matches_i(idx_ij_of_new_features_a,6)));
       count_1 = [];
       count_2 = [];
       vm = [];
       add_new_match = [];
       if ~isempty(missMatch)
           for m = 1:length(missMatch)
               if sum(PoseGraphMatrix(idx_ij_of_new_features_b(missMatch(m)),:)~=0) >= 3
                   % the index in new image needed to be delete 
                   count_1 = [count_1;missMatch(m)];
               else
                   % the point in all image needed to be delete
                   count_2 = [count_2;missMatch(m)];
                   vm = [vm;idx_ij_of_new_features_a(missMatch(m))];
               end
           end
           if ~isempty(count_1)
               PoseGraphMatrix(idx_ij_of_new_features_b(count_1),img_idx) = 0;
           end
           if ~isempty(count_2)
               % Since we need to eliminate some entire rows, we need to
               % delete the numbers in pointerToImages.
%                pointer_reduce = pointer_reduce + length(count_2);
               PoseGraphMatrix(idx_ij_of_new_features_b(count_2),:) = 0;
           end
           
           % record the 2D feature idx in count_2 for match_image, add it with current image as new feature in the following step
           add_new_match{i} = [matches_i(idx_ij_of_new_features_a(count_2),3) matches_i(idx_ij_of_new_features_a(count_2),6)];
           
           % clear 
           idx_ij_of_new_features_b([count_1;count_2]) = [];
           idx_ij_of_new_features_a([count_1;count_2]) = [];
       end
       PoseGraphMatrix(idx_ij_of_new_features_b, img_idx) = matches_i(idx_ij_of_new_features_a,6);
    end
    
    % 3. Taking new features in current image that haven't seen before      
    idx_features_imagei_seen_b4_based_preimgs = PoseGraphMatrix(1:pointerToImages(img_idx),img_idx) ~= 0;
    features_imagei_seen_b4 = PoseGraphMatrix(idx_features_imagei_seen_b4_based_preimgs,img_idx);
    % Adding these new features into PoseGraphMatrix    
    new_features_Imagei = setdiff(1:size(featureExtracted{img_idx},2), features_imagei_seen_b4);
    PoseGraphMatrix(1+pointerToImages(img_idx):length(new_features_Imagei)+pointerToImages(img_idx),img_idx) = new_features_Imagei';
    
    % add the miss match features before as new features
    if ~isempty(add_new_match)
       for i = 1:length(add_new_match)
           if ~isempty(add_new_match{i})
              new_match_id = zeros(size(add_new_match{i},1),1);
              for bi = 1:size(add_new_match{i},1)
                  new_match = find(PoseGraphMatrix(:,img_idx)==add_new_match{i}(bi,2));
                  if length(new_match) > 1
                     num_of_match = sum(PoseGraphMatrix(new_match,:)~=0,2); % if there are more than one row in new_match, find the row which viewed by most frames
                     [~,I] = max(num_of_match);
                     new_match_id(bi) = (I);
                     PoseGraphMatrix(setdiff(new_match,new_match(I)),:) = 0; % set other rows to zero
                  else
                     new_match_id(bi) = new_match;
                  end
              end
              PoseGraphMatrix(new_match_id,match_image_id(i)) = add_new_match{i}(:,1);
           end
        end
    end
    
end



        % 4. find features in image i that have been seen before 
        % clean up and look for agreement, otherwise it's a potential outlier and
        % should be removed!
        i = img_idx;
        if i > 1                    
            for l = unique(features_imagei_seen_b4)'
                % find the 3D index of the features seen before
                vp = find(PoseGraphMatrix(:,i) == l);
                % unify this common features that are seen in different
                % places but not pairwise common
                if(length(vp) > 1)
                % 1. Find the one which seen by most images
                    [~, most_poses] = max(sum(PoseGraphMatrix(vp,:)~=0,2));
                % 2. Extract the image index which can view this feature
                % except current image
                    nz_indices = find(PoseGraphMatrix(vp(most_poses),:) ~= 0);
                    nz_indices = setdiff(nz_indices, i);
                % 3. Scan through rest, if disagree, remove, else unify
                    % extract the 3D index of rest features (except most_poses)
                    potential_unifying_rows = setdiff(1:length(vp),most_poses);
                    % extract 2D features of potential_unifying_rows for
                    % each image
                    coeff_matrix = PoseGraphMatrix(vp(potential_unifying_rows), nz_indices);
                    % Either 0
                    status_matrix1 = (coeff_matrix == 0);
                    % Or equal to most_pose row
                    status_matrix2 = (coeff_matrix == repmat(PoseGraphMatrix(vp(most_poses),nz_indices), [length(vp)-1,1]));
                    % find the element which is neither 0 nor the most_pose
                    status_matrix = ~(status_matrix1 | status_matrix2);                    
                    
                    
                    % loop PoseGraphMatrix for the feature: 1) it is seen by current image; 2) it is either
                    % 0 or the most_pose by rest preview images
                    for kkk = potential_unifying_rows(sum(status_matrix,2) == 0)
                        % find the image index which can view the feature
                        % kkk except current image
                        nz_indices_kkk = setdiff(find(PoseGraphMatrix(vp(kkk), :) ~= 0),i);
                        for kkkk = nz_indices_kkk
                            ss_row1 = (PoseGraphMatrix(vp(setdiff(1:length(vp),most_poses)),kkkk) == 0);
                            ss_row2 = (PoseGraphMatrix(vp(setdiff(1:length(vp),most_poses)),kkkk) == PoseGraphMatrix(vp(kkk),kkkk));
                            ss_row = ~(ss_row1 | ss_row2);
                            if sum(ss_row) == 0
                               PoseGraphMatrix(vp(most_poses), kkkk) =  PoseGraphMatrix(vp(kkk), kkkk);
                            end
                        end
                    end
                    % Clear up after unifying
                    PoseGraphMatrix(vp(setdiff(1:length(vp),most_poses)),:) = 0;
                end
            end
        end

%     % clean up all zero rows
    matches_row = find( ...
        sum(PoseGraphMatrix(pointerToImages(i)+1:pointerToImages(i) + length(new_features_Imagei), i), 2) ~= 0);
    
    % extract only rows with available matches
    PoseGraphMatrix(1+pointerToImages(i):length(matches_row)+pointerToImages(i), :) = ...
        PoseGraphMatrix(matches_row + pointerToImages(i), :);
    
    % clean the features viewed by only one image
    mess_b4 = [];
    if  flag == 1
        mess_b4_candid = find(sum(PoseGraphMatrix(:,1:img_idx/2)~=0,2) == 1);
        for m = 1:length(mess_b4_candid)
            img_biew_m = find(PoseGraphMatrix(mess_b4_candid(m),:)); % check if this measurements are also single in the whole PoseGraphMatrix (they may have common feature in img_idx/2+1:img_idx)
            if length(img_biew_m) == 1
               mess_b4 = [mess_b4;mess_b4_candid(m)];
            end
        end
        PoseGraphMatrix(mess_b4,:) = 0;
    end
        
    
    
    % Update pointers
    pointer_clean = find(sum(PoseGraphMatrix(1:pointerToImages(i),:),2)==0);
    feature_view_b4_delete = pointer_clean;
    pointerToImages(i) = pointerToImages(i) - length(pointer_clean);
    pointerToImages(i+1) =  pointerToImages(i) + length(matches_row);

    

% find the index of features which viewed before but need to be deleted
% now. These features also need to be delete in FeaturesBag
% clear rows which are all zeros
mess_row = sum(PoseGraphMatrix,2)==0;
PoseGraphMatrix(mess_row,:) = [];
PoseGraphMatrix(pointerToImages(img_idx+1)+1:end, :) = []; 