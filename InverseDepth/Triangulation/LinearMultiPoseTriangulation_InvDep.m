function [temp_feat,temp_feat_xyz] = LinearMultiPoseTriangulation_InvDep(feat_id, ...
                                                 AbsolutePoses,         ...
                                                 allPoses, ...
                                                 PoseGraphMatrix,       ...
                                                 featureExtracted, ...                                         
                                                 CameraParams,CameraModel)
    
                                             
    posesViewkfeatIdx = find(PoseGraphMatrix(feat_id,allPoses));    
    
    if length(posesViewkfeatIdx) < 2 && size(PoseGraphMatrix,2) > 2
        temp_feat = [];
        temp_feat_xyz = [];
        return;
    end
    
    if length(posesViewkfeatIdx) < 3 && size(PoseGraphMatrix,2) > 5
        temp_feat = [];
        temp_feat_xyz = [];
        return;
    end                                         
    count = 0;
    A = zeros(3*length(posesViewkfeatIdx), 3);
    r = zeros(3*length(posesViewkfeatIdx), 1);
    % respect to the first camera view current feature
    ref = min(posesViewkfeatIdx);
      
    for l = posesViewkfeatIdx
        featureIdxImagel = PoseGraphMatrix(feat_id, allPoses(l));
        Cl_T_Cn = AbsolutePoses(:,:,allPoses(l))*[AbsolutePoses(:,1:3,ref)' -AbsolutePoses(:,1:3,ref)'*AbsolutePoses(:,4,ref);zeros(1,3) 1];
        bi = [StaticUndistort(CameraParams.fc,CameraParams.cc,CameraParams.kc,featureExtracted{allPoses(l)}(1:2,featureIdxImagel),CameraModel);1];
        bi = bi./norm(bi);
        
        A(3*count+1:3*count+3, :) = skewsymm(bi)*Cl_T_Cn(:,1:3);
        r(3*count+1:3*count+3, :) = -skewsymm(bi)*Cl_T_Cn(:,4);    
        count = count + 1;
    end
    ATA = A'*A;
    c = cond(ATA,2); % check the condition number of A
    
    if c >= 500 %2e8
        temp_feat = [];
        temp_feat_xyz = [];
        return;
    end
    V = ATA\(A'*r);
    
    rho = 1/norm(V(1:3));
    theta = atan2(V(2),V(1));
    phi = atan2(V(3),(V(1)^2+V(2)^2)^0.5);

    if rho > 5 || rho < 1e-2
        temp_feat= [];
        temp_feat_xyz = [];
        return;
    end
    
    temp_feat = [theta;phi;rho;ref];

    temp_feat_xyz = zeros(4,1);
    temp_feat_xyz(1:3) =  1/temp_feat(3)*[cos(temp_feat(2))*cos(temp_feat(1));cos(temp_feat(2))*sin(temp_feat(1));sin(temp_feat(2))];
    temp_feat_xyz(4) =  temp_feat(4);
    
end