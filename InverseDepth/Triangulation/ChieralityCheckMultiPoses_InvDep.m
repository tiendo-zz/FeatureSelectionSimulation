function infrontof = ChieralityCheckMultiPoses_InvDep(feat, AbsolutePoses, poses)

infrontof = 1;

if feat(3) > 5
    infrontof = 0;
    return;
end

feat_xyz = zeros(4,1);
feat_xyz(1:3) =  1/feat(3)*[cos(feat(2))*cos(feat(1));cos(feat(2))*sin(feat(1));sin(feat(2))];
feat_xyz(4) =  feat(4);

for k = poses
    Rti = AbsolutePoses(:,:,k)*[AbsolutePoses(:,1:3,feat_xyz(4))' -AbsolutePoses(:,1:3,feat_xyz(4))'*AbsolutePoses(:,4,feat_xyz(4));zeros(1,3) 1];
%     Rti = AbsolutePoses(:,:,k);
    Ci_feat = Rti*[feat_xyz(1:3);1];
    if(Ci_feat(3) <= 1e-3)
        infrontof = 0;
        break;
    end
end
