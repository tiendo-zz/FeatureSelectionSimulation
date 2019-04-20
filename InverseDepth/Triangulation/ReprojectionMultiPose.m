function triangulate_decision = ReprojectionMultiPose(feat_id, ...
                                                      PoseGraphMatrix,...
                                                      allPoses,...
                                                      AbsolutePoses,...
                                                      featureExtracted,...
                                                      temp_feat_xyz,...
                                                      CameraParams,...]
                                                      CameraModel,threshold)
                                                  
        posesViewkfeatIdx = find(PoseGraphMatrix(feat_id,allPoses)); 
        ref = min(posesViewkfeatIdx);
        % Check reprojection error:
        triangulate_decision = 1;
        for l = posesViewkfeatIdx
            featureIdxImagel = PoseGraphMatrix(feat_id, allPoses(l));
            Cl_T_Cn = AbsolutePoses(:,:,allPoses(l))*[AbsolutePoses(:,1:3,ref)' -AbsolutePoses(:,1:3,ref)'*AbsolutePoses(:,4,ref);zeros(1,3) 1];
            Cl_p_f = Cl_T_Cn*[temp_feat_xyz(1:3);1];
            Cl_p_f = Cl_p_f(1:3)/Cl_p_f(3);
            Cl_p_f = StaticDistort(CameraParams.fc,CameraParams.cc,CameraParams.kc,Cl_p_f(1:2),CameraModel);

            reproj_error = norm(Cl_p_f(1:2) - ...
                featureExtracted{allPoses(l)}(1:2,featureIdxImagel));
            if reproj_error > threshold
                triangulate_decision = 0;
                return;
            end
        end
end