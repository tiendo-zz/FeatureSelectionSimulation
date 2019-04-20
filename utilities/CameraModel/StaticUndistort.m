function homo = StaticUndistort(fc,cc,kc,z,CameraModel)

% change the pixel value to homo value
if isempty(kc)
    Kinv = [1/fc(1) 0 -cc(1)/fc(1); ...
        0 1/fc(2) -cc(2)/fc(2); ...
        0 0 1];
    homo = Kinv(1:2,:)*[z;ones(1,size(z,2))];
   
    return;
end

if strcmp(CameraModel,'Tango')
   homo = StaticUndistortTango(fc,cc,kc,z);
elseif strcmp(CameraModel,'Radial')
   homo = StaticUndistortRadial(fc,cc,kc,z);
elseif strcmp(CameraModel,'Fisheye')
   homo = StaticUndistortFisheye(fc,cc,kc,z);
end



