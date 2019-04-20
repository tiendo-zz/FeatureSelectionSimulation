function [AbsolutePoses,FeaturesBag,measurePoses,CameraParams,error_intrinsic] = simulation_full_BA_InvDep_pixel_2(PoseGraphMatrix,AbsolutePoses,featureExtracted,NumOfPoses,CameraParams,K_true)

error_intrinsic = zeros(4,NumOfPoses);
for newPose = 1:NumOfPoses
    if newPose==1
       AbsolutePoses(:,:,1) = [eye(3) zeros(3,1)];
       measurePoses = 1;
    else
       if newPose==2
         %% 5pt_ransac
         commen_feature_c1_c2 = find(PoseGraphMatrix(:,1)&PoseGraphMatrix(:,2));
         PairWiseMatches = zeros(length(commen_feature_c1_c2),6);
         PairWiseMatches(:,3) = PoseGraphMatrix(commen_feature_c1_c2,1);
         PairWiseMatches(:,6) = PoseGraphMatrix(commen_feature_c1_c2,2);
         PairWiseMatches(:,1:2) = featureExtracted{1}(1:2,PairWiseMatches(:,3))';
         PairWiseMatches(:,4:5) = featureExtracted{2}(1:2,PairWiseMatches(:,6))';
         PairWiseMatches_homo = PairWiseMatches;
         homo_1 = CameraParams.Kinv*[PairWiseMatches(:,1:2)';ones(1,size(PairWiseMatches(:,1:2),1))];
         PairWiseMatches_homo(:,1:2) = homo_1(1:2,:)';
         homo_2 = CameraParams.Kinv*[PairWiseMatches(:,4:5)';ones(1,size(PairWiseMatches(:,4:5),1))];
         PairWiseMatches_homo(:,4:5) = homo_2(1:2,:)';
         
         C2_R_C1 = AbsolutePoses(:,1:3,2);
         C2_t_C1 = AbsolutePoses(:,4,2);
         % Triangulation
%          FeaturesBag = LinearTriangulation_5pt_InvDep(PairWiseMatches_homo,PoseGraphMatrix, C2_R_C1', -C2_R_C1'*C2_t_C1);
         refPose = 1;
         scalePose = 2;
         else
           % change pixel value to homo 
           [Cn_R_Cr, Cn_t_Cr, InlierIdx, ~, Cr_P, Cn_z,idx_3d] = simulation_P3P_RANSAC_InvDep_pixel(PoseGraphMatrix, AbsolutePoses,  ...
                                                         featureExtracted,                      ...
                                                         measurePoses, newPose,                      ...
                                                         FeaturesBag, CameraParams,          ...
                                                         500, 0.7, 3);
           
           fprintf('P3P corresponding : %f\n', length(idx_3d)); 
           if  length(InlierIdx) < 15 
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
          [Cn_R_Cr, Cn_t_Cr, CameraParams, normHist] = simulation_PnP_NL_InvDep_pixel(Cn_R_Cr, Cn_t_Cr,AbsolutePoses, Cr_P(:, InlierIdx), Cn_z(InlierIdx,:)',CameraParams);                  
          AbsolutePoses(:,:,newPose) = [Cn_R_Cr, Cn_t_Cr];
          fprintf('Inlier new pose I%d: %d\n', newPose, length(InlierIdx));   

       end
       
      %% 7. Triangulation
      FeaturesBag = zeros(4, size(PoseGraphMatrix,1));
      [FeaturesBag] = simulation_MultiPoseTriangulation_InvDep_pixel(FeaturesBag,           ...
                                         AbsolutePoses,         ...
                                         PoseGraphMatrix,       ...
                                         newPose, measurePoses,    ...
                                         featureExtracted,      ...                                         
                                         CameraParams,3);

          
      fprintf('Number of triangulated features: %d\n', length(find(FeaturesBag(4,:) ~= 0)));
          
      %% 8. Bundle Adjustment
      measurePoses = [measurePoses newPose];
          
       % there need at least 3 images for ba   
       if length(measurePoses) > 2
%             [FeaturesBag, AbsolutePoses,CameraParams, ~, ~,~,~] = simulation_BA_InverseDepth_2viewBLS_Robust_pixel(FeaturesBag, AbsolutePoses, PoseGraphMatrix, featureExtracted, measurePoses, CameraParams,1e-2);
          [FeaturesBag, AbsolutePoses, CameraParams, ~, ~,~,~] = BA_InverseDepth_2viewBLS_pixel_simulation(FeaturesBag, AbsolutePoses, PoseGraphMatrix, featureExtracted, measurePoses, CameraParams);
%           [FeaturesBag, AbsolutePoses, ~, ~,~,~] = BA_InverseDepth_simulation(FeaturesBag, AbsolutePoses, PoseGraphMatrix, featureExtracted, measurePoses, CameraParams);
       end       
    end
    
    if newPose > 1
    for i = 1:newPose
        subplot(newPose,1,i); hold on
        plot(featureExtracted{i}(1,:),featureExtracted{i}(2,:),'r*');hold on
%         axis([0 CameraParams.cc(1)*2 0 CameraParams.cc(2)*2]);hold on
        triangulate_feat = find(FeaturesBag(4,:) & PoseGraphMatrix(:,i)');
        repro = zeros(3,length(triangulate_feat));
        for k = triangulate_feat
            m_theta_phi = [cos(FeaturesBag(2,k))*cos(FeaturesBag(1,k));cos(FeaturesBag(2,k))*sin(FeaturesBag(1,k));sin(FeaturesBag(2,k))];
            Ck_P_W = AbsolutePoses(:,:,i)*[InversePose(AbsolutePoses(:,:,FeaturesBag(4,k)));0 0 0 1]*[m_theta_phi;FeaturesBag(3,k)];
%             Ck_P_W = AbsolutePoses(:,:,i)*[InversePose(AbsolutePoses(:,:,FeaturesBag(4,k)));0 0 0 1]*[m_theta_phi;1];
            Ck_P_W = Ck_P_W/Ck_P_W(3);
            repro(:,k) = CameraParams.K*Ck_P_W;
            plot(repro(1,k),repro(2,k),'bo');hold on
        end
        
    end
    end
    fprintf('Image: %d\n', newPose);
    close all; 
    
    
    error_intrinsic(1,newPose) = abs(CameraParams.fc(1) - K_true(1,1));
    error_intrinsic(2,newPose) = abs(CameraParams.fc(2) - K_true(2,2));
    error_intrinsic(3,newPose) = abs(CameraParams.cc(1) - K_true(1,3));
    error_intrinsic(4,newPose) = abs(CameraParams.cc(2) - K_true(2,3));
end