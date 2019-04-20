function [R_true,t_true,idx]=TruePose_P3P_pixel(R,t,p,b,CameraParams,CameraModel)

%%%%%%%% input %%%%%%%%%%%%%%%%%%%%%%
% p: 3D features, w.r.t global frame
% b: 2D features, 3X1 vector, w.r.t camera frame
% K: intrinsic matrix
% R,t : w.r.t camera
% R_true, t_true : w.r.t camera

if size(b,1) ~= 3
    b = [b;ones(1,size(b,2))];
end

error = zeros(size(R,3),1);
% b = b/norm(b);
for i=1:size(R,3)
    homo = R(:,:,i)*p + t(:,i);
    homo = homo./homo(3);
    b_est = [StaticDistort(CameraParams.fc,CameraParams.cc,CameraParams.kc,homo,CameraModel);1];
    error(i) = norm(b_est-b);
end
[~,idx] = min(error);
R_true = R(:,:,idx);
t_true = t(:,idx);