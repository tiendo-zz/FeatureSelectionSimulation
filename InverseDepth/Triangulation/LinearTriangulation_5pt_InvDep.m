function Cr_P = LinearTriangulation_5pt_InvDep(matches, matches_pixel, PoseGraphMatrix, C1_R_C2, C1_t_C2,threshold,CameraParams)

global CameraModel

Y = LinearTriangulation_InverseDepth(matches, [eye(3) zeros(3,1)], [C1_R_C2(:,:,1)', -C1_R_C2(:,:,1)'*C1_t_C2]);

% Y_2 = LinearTriangulation(matches,[eye(3) zeros(3,1)],[C1_R_C2(:,:,1)', -C1_R_C2(:,:,1)'*C1_t_C2]);

Y_xyz = zeros(size(Y));
for k = 1:size(Y,2)
    Y_xyz(:,k) =  1/Y(3,k)*[cos(Y(2,k))*cos(Y(1,k));cos(Y(2,k))*sin(Y(1,k));sin(Y(2,k))];
end

ind = ChieralityCheck(Y_xyz, eye(3), zeros(3,1), C1_R_C2(:,:,1), C1_t_C2);
% ond = setdiff(1:size(Y,2), ind);


for l = ind
    C2_T_C1 = [C1_R_C2(:,:,1)', -C1_R_C2(:,:,1)'*C1_t_C2];
    C1_p_f = Y_xyz(1:3,l);
    C1_p_f = C1_p_f(1:3)/C1_p_f(3);
    C1_p_f = StaticDistort(CameraParams.fc,CameraParams.cc,CameraParams.kc,C1_p_f(1:2),CameraModel);
    C2_p_f = C2_T_C1*[Y_xyz(1:3,l);1];
    C2_p_f = C2_p_f(1:3)/C2_p_f(3);
    C2_p_f = StaticDistort(CameraParams.fc,CameraParams.cc,CameraParams.kc,C2_p_f(1:2),CameraModel);
    
    reproj_error_1 = norm(C1_p_f(1:2) - matches_pixel(l,1:2)');
    reproj_error_2 = norm(C2_p_f(1:2) - matches_pixel(l,4:5)');
    if reproj_error_1 > threshold || reproj_error_2 > threshold
       ind(find(ind==l)) = 0;
    end
end

ind(find(ind==0)) = [];

Cr_P = zeros(4, size(PoseGraphMatrix,1));
[update_indices_a, update_indices] = ismember(matches(ind,3)', PoseGraphMatrix(:,1));
update_indices_a = find(update_indices_a~=0);
update_indices = update_indices(update_indices~=0);
Cr_P(:,update_indices) = [Y(:, ind(update_indices_a)); 1*ones(1,length(update_indices_a))];
