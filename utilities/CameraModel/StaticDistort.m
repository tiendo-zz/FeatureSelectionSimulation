function  z = StaticDistort(fc,cc,kc,homo,CameraModel)


% change the homo value to pixel value
if isempty(kc)
    K = [fx 0 cx; ...
        0 fy cy; ...
        0 0 1];
    z = K(1:2,:)*[homo;ones(1,size(homo,2))];

    return;
end

if strcmp(CameraModel,'Tango')
   z = StaticDistortTango(fc,cc,kc,homo);
elseif strcmp(CameraModel,'Radial')
   z = StaticDistortRadial(fc,cc,kc,homo);
elseif strcmp(CameraModel,'Fisheye')
   z = StaticDistortFisheye(fc,cc,kc,homo);
end
