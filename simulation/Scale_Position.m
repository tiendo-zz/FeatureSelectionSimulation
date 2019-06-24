function [AbsolutePoses_InvDep,FeaturesBag_InvDep_xyz,scale_InvDep] = Scale_Position(AbsolutePoses_true,AbsolutePoses_InvDep,FeaturesBag_InvDep_xyz,FeatureBag_true_xyz)

scale_InvDep = norm(AbsolutePoses_true(:,4,2));
for i = 2:size(AbsolutePoses_InvDep,3)
    AbsolutePoses_InvDep(:,4,i) = AbsolutePoses_InvDep(:,4,i) * scale_InvDep;
end

for i = 1:size(FeaturesBag_InvDep_xyz, 2)
    FeaturesBag_InvDep_xyz(1:3, i) = FeaturesBag_InvDep_xyz(1:3, i) * scale_InvDep;
end