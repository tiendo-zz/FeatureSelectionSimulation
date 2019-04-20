function   [FeaturesBag] = simulation_MultiPoseTriangulation_InvDep_pixel(FeaturesBag,           ...
                                                 AbsolutePoses,         ...
                                                 PoseGraphMatrix,       ...
                                                 newPose, usedPoses,    ...
                                                 featureExtracted, ...                                         
                                                 CameraParams,threshold)
                                                                                         
                                             
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
    
    if length(posesViewkfeatIdx) < 2 && newPose > 2
        continue;
    end
    
    if length(posesViewkfeatIdx) < 3 && newPose > 5
        continue;
    end
    
    count = 0;
    A = zeros(3*length(posesViewkfeatIdx), 3);
    r = zeros(3*length(posesViewkfeatIdx), 1);
    % respect to the first camera view current feature
    ref = min(posesViewkfeatIdx);
    
%     if (FeaturesBag(4,k) ~= 0) && (ref ~= FeaturesBag(4,k))
%         fprintf('k: %d\n', k);
%     end
        
    
    
    for l = posesViewkfeatIdx
        featureIdxImagel = PoseGraphMatrix(k, allPoses(l));
        Cl_T_Cn = AbsolutePoses(:,:,allPoses(l))*[AbsolutePoses(:,1:3,ref)' -AbsolutePoses(:,1:3,ref)'*AbsolutePoses(:,4,ref);zeros(1,3) 1];
        bi = CameraParams.Kinv*[featureExtracted{allPoses(l)}(1:2,featureIdxImagel);1];
        bi = bi./norm(bi);
        
        A(3*count+1:3*count+3, :) = skewsymm(bi)*Cl_T_Cn(:,1:3);
        r(3*count+1:3*count+3, :) = -skewsymm(bi)*Cl_T_Cn(:,4);    
        count = count + 1;
    end
    ATA = A'*A;
    c = cond(ATA,2); % check the condition number of A
    
    if c >= 2e8
        continue;
    end
    V = ATA\(A'*r);
    
    rho = 1/norm(V(1:3));
    theta = atan2(V(2),V(1));
    phi = atan2(V(3),(V(1)^2+V(2)^2)^0.5);

    if rho > 5 || rho < 1e-2
        continue;
    end
    
    temp_feat = [theta;phi;rho;ref];

    temp_feat_xyz = zeros(4,1);
    temp_feat_xyz(1:3) =  1/temp_feat(3)*[cos(temp_feat(2))*cos(temp_feat(1));cos(temp_feat(2))*sin(temp_feat(1));sin(temp_feat(2))];
    temp_feat_xyz(4) =  temp_feat(4);
    
    if( ChieralityCheckMultiPoses(temp_feat_xyz, AbsolutePoses, allPoses(posesViewkfeatIdx)) )
        % Check reprojection error:
        triangulate_decision = 1;
        for l = posesViewkfeatIdx
            featureIdxImagel = PoseGraphMatrix(k, allPoses(l));
            Cl_T_Cn = AbsolutePoses(:,:,allPoses(l))*[AbsolutePoses(:,1:3,ref)' -AbsolutePoses(:,1:3,ref)'*AbsolutePoses(:,4,ref);zeros(1,3) 1];
            Cl_p_f = Cl_T_Cn*[temp_feat_xyz(1:3);1];
            Cl_p_f = Cl_p_f(1:3)/Cl_p_f(3);
            Cl_p_f = CameraParams.K*[Cl_p_f(1:2);1];

            reproj_error = norm(Cl_p_f(1:2) - ...
                featureExtracted{allPoses(l)}(1:2,featureIdxImagel));
            if reproj_error > threshold
                triangulate_decision = 0;
                break;
            end
        end
        if triangulate_decision == 1
            FeaturesBag(:,k) = temp_feat;
%             FeaturesBag_xyz(:,k) = temp_feat_xyz;
%         else
%             FeaturesBag(:,k) = zeros(4,1);
        end
    end
end


