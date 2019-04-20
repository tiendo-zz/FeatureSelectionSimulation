function [inlier_idx,error_mat] = ComputeInlier3D2DReprojection_pixel(P, w, CameraParams, R, t, threshold,CameraModel)
% R,t: w.r.t camera

w_ = [R t]*[P(1:3,:);ones(1,size(P,2))];
w_(1,:) = w_(1,:)./w_(3,:);
w_(2,:) = w_(2,:)./w_(3,:);
w_(3,:) = w_(3,:)./w_(3,:);
w_ = StaticDistort(CameraParams.fc,CameraParams.cc,CameraParams.kc,w_(1:2,:),CameraModel);
dw = w_' - w(:,1:2);
error_mat = sqrt(dw(:,1).*dw(:,1) + dw(:,2).*dw(:,2));

inlier_idx = find(error_mat < threshold);

