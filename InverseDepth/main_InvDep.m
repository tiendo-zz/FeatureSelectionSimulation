clear all
close all
clc
addpath('./Jacobian/');
addpath('./PnP/');
addpath('./Triangulation/');
addpath('./BundleAdjustment/');
addpath('../utilities/Drawings/');
addpath('../utilities/robotics3D/');
addpath('../utilities/CameraModel/');
addpath('../utilities/PoseGraphMatrix/');
addpath('../utilities/image-selection/');
addpath('../utilities/P3P/');
addpath('../TwoViewReconstruction/');
addpath('../RobustCostFunction');

% rng('default');
% rng(1);

global datapath image_folder CameraModel device output_folder config_folder
datapath = '/usr/local/google/home/feiwu/Fei/PTAM/';
% image_folder = '/usr/local/google/home/feiwu/Documents/datasets/vicon_test/fisheye';
config_folder = '/usr/local/google/home/feiwu/Documents/datasets/testmm/';
image_folder = '/usr/local/google/home/feiwu/Documents/datasets/testmm/images';
output_folder = [datapath,'img_data'];
CameraModel = "Radial";
device = "pixel_2";

% import new image by FOR loop, loop all image numbers
PoseGraphMatrix = [];
pointerToImages = [];


if strcmpi(device,'vicon')
   % vicon test
   fc = [254.997, 254.933];
   cc = [326.69, 230.118];
   kc = [0.925175, 0.0, 0.0, 0.0, 0.0];
   NumOfPoses = 1220;
elseif strcmpi(device,'pixel')
    % pixel phone
    fc = [479.6077759071951,479.3840971418472];
    cc = [321.584578057932,241.2778075680262]; 
    kc = [0.01756713993446139 -0.05496311956671585 0 0 0.08038811792685442]; 
    NumOfPoses = 1560;
elseif strcmpi(device,'pixel_2')
    fc = [476.7507599247749,  476.5919885539533];
    cc = [319.0116340608212,  239.7613723951213];
    kc = [0.02636362390446432,  -0.07160736234062902, 0.0, 0.0, 0.09555537909279475];
    NumOfPoses = 1600;
end

fx = fc(1); fy = fc(2);
cx = cc(1); cy = cc(2);
K    = [fx 0 cx; ...
        0 fy cy; ...
        0 0 1];
Kinv = [1/fx 0 -cx/fx; ...
        0 1/fy -cy/fy; ...
        0 0 1]; 
CameraParams = struct;
CameraParams.fc = fc; 
CameraParams.cc = cc; CameraParams.kc = kc; 
CameraParams.K = K; CameraParams.Kinv = Kinv;

% parameters
multiPosesVisualization = 0;
reprojVisualization = 0;
VisualizeMacthing = 0;
CreatePoseGraphMatrix = 1;
ApplyVT = 1;
FeaturesBag = [];

newPose = 0;
begin_idx = newPose + 1;

cur_flag_id_begin = 5;
cur_flag_id = cur_flag_id_begin;  
rate_init = 98;
rate = rate_init;

ImageId= 0;
% while ImageId < NumOfPoses
inlier_ratio = zeros(30,1);
while newPose < 100 && ImageId < NumOfPoses
    newPose = newPose + 1;
    if newPose == 1
       %% 1. Extract features
       command = [datapath,'ExtractFREAK.sh ',image_folder,' ',config_folder,' ',output_folder,' ',num2str(80),' ',num2str(newPose),' ',num2str(1)];
       tic;
       system(command);
       t(newPose,1) = toc;
       ImageId = load([output_folder,'/idx_selected.txt']);
           
       if ApplyVT == 1
         dlmwrite('measurePosesId.txt',[]);
       end
       
       AbsolutePoses(:,:,1) = [eye(3) zeros(3,1)];
       measurePoses = 1;
       measurePosesId = ImageId(1);
       
       command = [datapath,'QueryVT.sh ',image_folder,' ',config_folder,' ',output_folder,' ',num2str(ImageId),' ',num2str(1)];
       system(command);
       
       new_FREAK = load([output_folder,'/features_clone/',num2str(ImageId),'.txt']);
       featureExtracted{1} = new_FREAK(2:end,1:2)';
       

       featureExtracted_homo{1}(1:2,:) = StaticUndistort(fc,cc,kc,featureExtracted{1}(1:2,:),CameraModel);

       %% 2. add into PoseGraphMatrix
       PairWiseMatches = cell(nchoosek(2,2),1);  
       match_image_id = 0;
       if CreatePoseGraphMatrix == 1
          [PoseGraphMatrix, pointerToImages] = ConstructPoseGraphMatrix(featureExtracted, PairWiseMatches, PoseGraphMatrix, pointerToImages, newPose,match_image_id,0);
       else
           PoseGraphMatrix = load('PoseGraphMatrix.mat','PoseGraphMatrix');
           PoseGraphMatrix = PoseGraphMatrix.PoseGraphMatrix;
       end
       I_rgb{newPose} = imread([image_folder,'/m',num2str(ImageId,'%07d'),'.pgm']);
   else
       if length(measurePoses) < 4
          match_id_candid = measurePoses;
       else
          match_id_candid = measurePoses((length(measurePoses)-3):end); 
       end
       
       if ApplyVT == 1
           if newPose == cur_flag_id
               if newPose == cur_flag_id_begin
                   cur_flag = 2; % initial vt
               else
                   cur_flag = 1; % build VT tree
               end
           
               % update cur_flag_id
               if cur_flag_id < 100
                   cur_flag_id = cur_flag_id * 2;
               else
                   cur_flag_id = cur_flag_id + 100;
               end
           else
               cur_flag = 0; 
           end
       else
           cur_flag = 3; % no VT
       end
       
       [PairWiseMatches,PairWiseMatches_homo,featureExtracted,match_image_id,ImageId,t] = ImageSelect_5dofRANSAC(newPose,featureExtracted,match_id_candid,CameraParams,measurePosesId,rate,NumOfPoses,t);
       if isempty(match_image_id)
          if rate == 100
             return;
          end
          rate = rate + 2;
          featureExtracted{newPose} = [];
          featureExtracted_homo{newPose} = [];
          featureDescriptor{newPose} = [];
          newPose = newPose - 1;    
          continue;
       end
       
       % load image
       I_rgb{newPose} = imread([image_folder,'/m',num2str(ImageId,'%07d'),'.pgm']);
       
       
       %% 3. add new image into PoseGraphMatrix
       if CreatePoseGraphMatrix == 1
          PoseGraphMatrix_pre = PoseGraphMatrix;
          pointerToImages_pre = pointerToImages;
          FeaturesBag_pre = FeaturesBag;
          [PoseGraphMatrix, pointerToImages,feature_view_b4_delete] = ConstructPoseGraphMatrix(featureExtracted, PairWiseMatches, PoseGraphMatrix, pointerToImages,newPose,match_image_id,cur_flag);
          if ~isempty(FeaturesBag) 
             if ~isempty(feature_view_b4_delete)          
                FeaturesBag(:,feature_view_b4_delete) = [];
             end
             FeaturesBag_old = FeaturesBag;
             FeaturesBag = zeros(4,size(PoseGraphMatrix,1));
             FeaturesBag(:,1:size(FeaturesBag_old,2)) = FeaturesBag_old;
          end
       else
           PoseGraphMatrix = load('PoseGraphMatrix.mat','PoseGraphMatrix');
           PoseGraphMatrix = PoseGraphMatrix.PoseGraphMatrix;
       end
       
%        % display image and matching points
       if VisualizeMacthing == 1
          idx1 = match_id_candid(1);
          idx2 = newPose;
          I1 = I_rgb{measurePoses(match_id_candid(1))};
          I2 = I_rgb{newPose};
          matches = [];
          matches(:,1) = PoseGraphMatrix((PoseGraphMatrix(:,idx1)~=0)&(PoseGraphMatrix(:,idx2)~=0),idx1);
          matches(:,2) = PoseGraphMatrix((PoseGraphMatrix(:,idx1)~=0)&(PoseGraphMatrix(:,idx2)~=0),idx2);
          VisualizeMatches(I1, I2, matches, featureExtracted{idx1},featureExtracted{idx2},'b' ,'g');
          VisualizeMatches_2(I1, I2, PairWiseMatches{SUTrilInd(idx2,idx1,idx2)}, 'b' ,'g');
       end
    
    %% 4. P3P or 5-point-algorithm
       CameraParams.K = eye(3);
       CameraParams.Kinv = eye(3);
       if newPose==2
           AbsolutePoses(:,:,newPose) = load([output_folder,'/match_data/CameraPose_',num2str(measurePosesId(newPose-1)),'_',num2str(ImageId),'.txt']);
%            [AbsolutePoses(:,:,newPose),PairWiseMatches_homo{SUTrilInd(2,1,2)},c1_f_hat] = Two_View_BLS(PairWiseMatches_homo{SUTrilInd(2,1,2)},AbsolutePoses(:,:,newPose));
           
%            FeaturesBag = zeros(4, size(PoseGraphMatrix,1));
%            [update_indices_a, update_indices] = ismember(PairWiseMatches_homo{SUTrilInd(2,1,2)}(:,3)', PoseGraphMatrix(:,1));
%            update_indices_a = find(update_indices_a~=0);
%            update_indices = update_indices(update_indices~=0);
%            FeaturesBag(:,update_indices) = [c1_f_hat(:, update_indices_a); 1*ones(1,length(update_indices_a))];      
%            Pose_scale = 0.0441;
           Pose_scale = 1;
           AbsolutePoses(:,4,newPose) = AbsolutePoses(:,4,newPose) * Pose_scale; % change scale
           C2_R_C1 = AbsolutePoses(:,1:3,newPose);
           C2_t_C1 = AbsolutePoses(:,4,newPose);
           
           if strcmpi(CameraModel,'Tango')
              featureExtracted_homo{2}(1:2,:) = StaticUndistort(fc,cc,kc,featureExtracted{2}(1:2,:));
           elseif strcmpi(CameraModel,'Radial')
              featureExtracted_homo{2}(1:2,:) = StaticUndistortRadial(fc,cc,kc,featureExtracted{2}(1:2,:));
           end
           
           
           % Triangulation
           FeaturesBag = LinearTriangulation_5pt_InvDep(PairWiseMatches_homo{SUTrilInd(newPose,1,newPose)},PoseGraphMatrix, CameraParams, C2_R_C1', -C2_R_C1'*C2_t_C1);
           refPose = 1;
           scalePose = 2;
           
       else
           % change pixel value to homo 

           if strcmpi(CameraModel,'Tango')
              featureExtracted_homo{newPose}(1:2,:) = StaticUndistort(fc,cc,kc,featureExtracted{newPose}(1:2,:));
           elseif strcmpi(CameraModel,'Radial')
              featureExtracted_homo{newPose}(1:2,:) = StaticUndistortRadial(fc,cc,kc,featureExtracted{newPose}(1:2,:));
           end
           

           [Cn_R_Cr, Cn_t_Cr, InlierIdx, ~, Cr_P, Cn_z,idx_3d] = P3P_RANSAC_InvDep(PoseGraphMatrix, AbsolutePoses,  ...
                                                         featureExtracted_homo,                      ...
                                                         measurePoses, newPose,                      ...
                                                         FeaturesBag, eye(3),          ...
                                                         500, 0.7, 0.01);

           fprintf('P3P corresponding : %f\n', length(idx_3d)); 
           if  length(InlierIdx) < 15
               if rate == 100
                   return;
               end
               rate = rate + 2;
               PoseGraphMatrix = PoseGraphMatrix_pre;
               pointerToImages = pointerToImages_pre;
               FeaturesBag = FeaturesBag_pre;

               featureExtracted{newPose} = [];
               featureExtracted_homo{newPose} = [];
               featureDescriptor{newPose} = [];
               newPose = newPose - 1;
               
               continue;
           else
               inlier_ratio(newPose) = length(InlierIdx)/length(idx_3d);
               rate = rate_init;
           end
           
          
          % Visualization
          if reprojVisualization == 1     
             reprojection_matches = ReprojectionToMatches(Cr_P, Cn_R_Cr, Cn_t_Cr, CameraParams);
             ReprojectionVisualization2(I_rgb{newPose}, Cn_z, reprojection_matches); 
          end                                           
          
          
          
          %% 5. update pose graph matrix based p3p inliers
          [PoseGraphMatrix, pointerToImages,feature_view_b4_delete] = ConstructPoseGraphMatrix_p3p(InlierIdx, idx_3d, PoseGraphMatrix, pointerToImages, newPose, FeaturesBag);
          if ~isempty(FeaturesBag) 
             if ~isempty(feature_view_b4_delete)
                FeaturesBag(:,feature_view_b4_delete) = [];
             end
             FeaturesBag_old = FeaturesBag;
             FeaturesBag = zeros(4,size(PoseGraphMatrix,1));
             FeaturesBag(:,1:size(FeaturesBag_old,2)) = FeaturesBag_old;
          end
          
          
          %% 6. nonlinear PnP
          [Cn_R_Cr, Cn_t_Cr, normHist] = PnP_NL_InvDep(Cn_R_Cr, Cn_t_Cr,AbsolutePoses, Cr_P(:, InlierIdx), Cn_z(InlierIdx,:)',CameraParams);                      

          
          AbsolutePoses(:,:,newPose) = [Cn_R_Cr, Cn_t_Cr];
          fprintf('Inlier new pose I%d: %d\n', newPose, length(InlierIdx));   
          
          % Visualization
          if reprojVisualization == 1     
             reprojection_matches = ReprojectionToMatches(Cr_P, Cn_R_Cr, Cn_t_Cr, CameraParams);
             ReprojectionVisualization2(I_rgb{newPose}, Cn_z, reprojection_matches); 
          end  
       end
       
      %% 7. Triangulation
      FeaturesBag = MultiPoseTriangulation_InvDep(FeaturesBag,           ...
                                         AbsolutePoses,         ...
                                         PoseGraphMatrix,       ...
                                         newPose, measurePoses,    ...
                                         featureExtracted_homo,      ...                                         
                                         CameraParams,0.01);

          
      fprintf('Number of triangulated features: %d\n', length(find(FeaturesBag(4,:) ~= 0)));
          
      % Visualization
      if multiPosesVisualization == 1
         VisualizeMultiPoses(AbsolutePoses, FeaturesBag, measurePoses, newPose);
         VisualizeAllPosesReprojection(PoseGraphMatrix, AbsolutePoses, FeaturesBag, featureExtracted, measurePoses, newPose, I_rgb, CameraParams);
      end
          
      %% 8. Bundle Adjustment
      measurePoses = [measurePoses newPose];
      measurePosesId = [measurePosesId;ImageId];
          
       % there need at least 3 images for ba   
       if length(measurePoses) > 2
%           [FeaturesBag, AbsolutePoses, ~, ~,EraseId] = BA_InverseDepth_2viewBLS(FeaturesBag, AbsolutePoses, PoseGraphMatrix, featureExtracted_homo, measurePoses, CameraParams);
          [FeaturesBag, AbsolutePoses, ~, ~,~,~] = BA_InverseDepth_2viewBLS_Robust(FeaturesBag, AbsolutePoses, PoseGraphMatrix, featureExtracted_homo, measurePoses, CameraParams,1e-1);
       end
       
%        % record reprojection error of camera newPose
%        repro_pose = newPose;
%        feat_view_3d_idx = find(FeaturesBag(4,:) == repro_pose);
%        feat_view_2d_idx = PoseGraphMatrix(feat_view_3d_idx,repro_pose);
%        feat_3d = [FeaturesBag(1:3,feat_view_3d_idx);ones(1,length(feat_view_3d_idx))];
%        feat_2d = featureExtracted_homo{repro_pose}(:,feat_view_2d_idx);
%        feat_repro = feat_3d(1:3,:);
%        feat_repro(1:2,:) = feat_repro(1:2,:)./feat_repro(3,:);
%        repro_error = [repro_error;norm(feat_repro(1:2,:) - feat_2d)];
    

                   
    % add image into vt tree
    command = [datapath,'QueryVT.sh ',image_folder,' ',config_folder,' ',output_folder,' ',num2str(ImageId),' ',num2str(2)];
    system(command);

    end
    if ApplyVT == 1
       dlmwrite('measurePosesId.txt',measurePosesId);
    end
    
    
    fprintf('Image: %d\n', ImageId);
    close all;
    
    % save tmp dataset
    
%     save('tmp/AbsolutePoses.mat','AbsolutePoses');
%     save('tmp/FeaturesBag.mat','FeaturesBag');
%     save('tmp/PoseGraphMatrix.mat','PoseGraphMatrix');
%     save('tmp/featureExtracted.mat','featureExtracted');
%     save('tmp/pointerToImages.mat','pointerToImages');
%     save('tmp/featureExtracted_homo.mat','featureExtracted_homo');
%     save('tmp/measurePoses.mat','measurePoses');
%     save('tmp/measurePosesId.mat','measurePosesId'); 
    
       
end
% save('dataset/inlier_ratio_1.mat','inlier_ratio');


% remove outliers
featureAvailableList = find(FeaturesBag(4,:));
featureAvailableList = featureAvailableList-1;


% transform FeaturesBag of Inverse Depth param into Cartesian param
FeaturesBag_InvDep_xyz = zeros(size(FeaturesBag));
for k = 1:size(FeaturesBag,2)
    if FeaturesBag(4,k)~=0
       FeaturesBag_InvDep_xyz(1:3,k) =  1/FeaturesBag(3,k)*[cos(FeaturesBag(2,k))*cos(FeaturesBag(1,k));cos(FeaturesBag(2,k))*sin(FeaturesBag(1,k));sin(FeaturesBag(2,k))];
       FeaturesBag_InvDep_xyz(4,k) = FeaturesBag(4,k);
    end
end



% save data
if ApplyVT == 0
       save('dataset/AbsolutePoses_NoVT.mat','AbsolutePoses');
       save('dataset/FeaturesBag_NoVT.mat','FeaturesBag');
       save('dataset/PoseGraphMatrix_NoVT.mat','PoseGraphMatrix');
       save('dataset/featureExtracted_NoVT.mat','featureExtracted');
       save('dataset/pointerToImages_NoVT.mat','pointerToImages');
       save('dataset/featureExtracted_homo_NoVT.mat','featureExtracted_homo');
       save('dataset/FeaturesBag_InvDep_xyz_NoVT.mat','FeaturesBag_InvDep_xyz');
       save('dataset/measurePoses_NoVT.mat','measurePoses');
       save('dataset/measurePosesId_NoVT.mat','measurePosesId'); 
else
       save('dataset/AbsolutePoses_VT.mat','AbsolutePoses');
       save('dataset/FeaturesBag_VT.mat','FeaturesBag');
       save('dataset/PoseGraphMatrix_VT.mat','PoseGraphMatrix');
       save('dataset/featureExtracted_VT.mat','featureExtracted');
       save('dataset/pointerToImages_VT.mat','pointerToImages');
       save('dataset/featureExtracted_homo_VT.mat','featureExtracted_homo');
       save('dataset/FeaturesBag_InvDep_xyz_VT.mat','FeaturesBag_InvDep_xyz');
       save('dataset/measurePoses_VT.mat','measurePoses');
       save('dataset/measurePosesId_VT.mat','measurePosesId'); 
end


% display
% figure;
% plot(repro_error);
% title('reprojection error');
% if ApplyVT == 1
%      saveas(gcf,[datapath,'Output/repro_error_VT.jpg']);
% else
%     saveas(gcf,[datapath,'Output/repro_error_NoVT.jpg']);
% end

figure;
for i = 1:newPose
plot3(AbsolutePoses(1,4,i),AbsolutePoses(2,4,i),AbsolutePoses(3,4,i),'r*');hold on
end
title('camera pose');
if ApplyVT == 1
     saveas(gcf,[datapath,'Output/CameraPosition_VT.fig']);
else
    saveas(gcf,[datapath,'Output/CameraPosition_NoVT.fig']);
end

% figure;
% plot(t(:,2));hold on
% plot(t(:,3));
% legend('ba','pnp');
% title('time cost');
% saveas(gcf,[datapath,'Output/time cost.jpg']);

VisualizeMultiPoses(AbsolutePoses, FeaturesBag_InvDep_xyz, measurePoses(1:end-1), newPose);
if ApplyVT == 1
    saveas(gcf,[datapath,'Output/VT_MuliPoses.fig']);
else
    saveas(gcf,[datapath,'Output/NoVT_MuliPoses.fig']);
end
% CameraParams.fc = fc; 
% CameraParams.cc = cc; CameraParams.kc = kc; 
% CameraParams.K = K; CameraParams.Kinv = Kinv;
% ShowPose = [1;10;20;30];
% measurePoses = measurePoses(ShowPose);
% measurePosesId = measurePosesId(ShowPose);
% VisualizeAllPosesReprojection(PoseGraphMatrix, AbsolutePoses, FeaturesBag_InvDep_xyz, featureExtracted, measurePoses(1:end-1), measurePoses(end), measurePosesId, I_rgb, CameraParams);
% if ApplyVT == 1
%     saveas(gcf,[datapath,'Output/VT_AllPosesReprojection.jpg']);
% else
%     saveas(gcf,[datapath,'Output/NoVT_AllPosesReprojection.jpg']);
% end
% 