function [inlier_idx] = simulation_ComputeInlier3D2DReprojection_pixel(P, w, CameraParams, R, t, threshold)
% R,t: w.r.t camera

w_ = [R t]*[P(1:3,:);ones(1,size(P,2))];
w_(1,:) = w_(1,:)./w_(3,:);
w_(2,:) = w_(2,:)./w_(3,:);
w_(3,:) = w_(3,:)./w_(3,:);
w_ = CameraParams.K*[w_(1:2,:);ones(1,size(w_,2))];
dw = w_(1:2,:)' - w(:,1:2);
error_mat = sqrt(dw(:,1).*dw(:,1) + dw(:,2).*dw(:,2));

inlier_idx = find(error_mat < threshold);

