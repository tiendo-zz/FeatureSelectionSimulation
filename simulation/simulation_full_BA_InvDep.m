function [AbsolutePoses,FeaturesBag,measurePoses] = simulation_full_BA_InvDep(PoseGraphMatrix,AbsolutePoses,featureExtracted,NumOfPoses)

for newPose = 1:NumOfPoses
    if newPose==1
       AbsolutePoses(:,:,1) = [eye(3) zeros(3,1)];
       measurePoses = 1;
    else
       CameraParams.K = eye(3);
       CameraParams.Kinv = eye(3);
       if newPose==2
         %% 5pt_ransac
         commen_feature_c1_c2 = find(PoseGraphMatrix(:,1)&PoseGraphMatrix(:,2));
         PairWiseMatches = zeros(length(commen_feature_c1_c2),6);
         PairWiseMatches(:,3) = PoseGraphMatrix(commen_feature_c1_c2,1);
         PairWiseMatches(:,6) = PoseGraphMatrix(commen_feature_c1_c2,2);
         PairWiseMatches(:,1:2) = featureExtracted{1}(1:2,PairWiseMatches(:,3))';
         PairWiseMatches(:,4:5) = featureExtracted{2}(1:2,PairWiseMatches(:,6))';
         
         % normalize 2nd pose translation
         AbsolutePoses(:,4,2) = AbsolutePoses(:,4,2)/norm(AbsolutePoses(:,4,2));
         C2_R_C1 = AbsolutePoses(:,1:3,2);
         C2_t_C1 = AbsolutePoses(:,4,2);
         % Triangulation

         FeaturesBag = LinearTriangulation_5pt_InvDep_nointrinsic(PairWiseMatches, ...
                                                    PoseGraphMatrix, ...
                                                    C2_R_C1', ...
                                                    -C2_R_C1'*C2_t_C1, 1e-1);

                                                  
         refPose = 1;
         scalePose = 2;
         else
           % change pixel value to homo 
           [Cn_R_Cr, Cn_t_Cr, InlierIdx, ~, Cr_P, Cn_z, idx_3d] = P3P_RANSAC_InvDep(PoseGraphMatrix, AbsolutePoses,  ...
                                                         featureExtracted,                      ...
                                                         measurePoses, newPose,                      ...
                                                         FeaturesBag, eye(3),          ...
                                                         1000, 0.7, 1e-3);
           
           fprintf('P3P corresponding : %f\n', length(idx_3d)); 
           if  length(InlierIdx) < 4 
               return;
           end
                                       
                   
          %% 5. update pose graph matrix based p3p inliers
%           [PoseGraphMatrix, pointerToImages,feature_view_b4_delete] = ConstructPoseGraphMatrix_p3p(InlierIdx, idx_3d, PoseGraphMatrix, pointerToImages, newPose);
%           if ~isempty(FeaturesBag) 
%              if ~isempty(feature_view_b4_delete)
%                 FeaturesBag(:,feature_view_b4_delete) = [];
%              end
%              FeaturesBag_old = FeaturesBag;
%              FeaturesBag = zeros(4,size(PoseGraphMatrix,1));
%              FeaturesBag(:,1:size(FeaturesBag_old,2)) = FeaturesBag_old;
%           end
          
          
          %% 6. nonlinear PnP
          [Cn_R_Cr, Cn_t_Cr, normHist] = PnP_NL_InvDep(Cn_R_Cr, Cn_t_Cr,AbsolutePoses, Cr_P(:, InlierIdx), Cn_z(InlierIdx,:)',CameraParams);                  
          AbsolutePoses(:,:,newPose) = [Cn_R_Cr, Cn_t_Cr];
          fprintf('Inlier new pose I%d: %d\n', newPose, length(InlierIdx));   

       end
       
      %% 7. Triangulation
      if length(measurePoses) > 2  
        FeaturesBag = MultiPoseTriangulation_InvDep_nointrinsic(FeaturesBag,           ...
                                                              AbsolutePoses,         ...
                                                              PoseGraphMatrix,       ...
                                                              newPose, measurePoses,    ...
                                                              featureExtracted,      ...                                         
                                                              CameraParams,0.01);
      end

          
      fprintf('Number of triangulated features: %d\n', length(find(FeaturesBag(4,:) ~= 0)));
          
      %% 8. Bundle Adjustment
      measurePoses = [measurePoses newPose];
          
       % there need at least 3 images for ba   
       if length(measurePoses) > 2           
           [FeaturesBag, AbsolutePoses, ~, ~,~] = BA_InverseDepth_2viewBLS(FeaturesBag, AbsolutePoses, PoseGraphMatrix, featureExtracted, measurePoses, CameraParams);
%           [FeaturesBag, AbsolutePoses, ~, ~,~,~] = BA_InverseDepth_2viewBLS_simulation(FeaturesBag, AbsolutePoses, PoseGraphMatrix, featureExtracted, measurePoses, CameraParams);
%           [FeaturesBag, AbsolutePoses, ~, ~,~,~] = BA_InverseDepth_simulation(FeaturesBag, AbsolutePoses, PoseGraphMatrix, featureExtracted, measurePoses, CameraParams);
       end       
    end
    
    fprintf('Image: %d\n', newPose);
    close all;      
end