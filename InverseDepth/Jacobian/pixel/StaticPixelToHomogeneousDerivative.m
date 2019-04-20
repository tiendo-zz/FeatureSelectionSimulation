function H_pixel = StaticPixelToHomogeneousDerivative(undistorted_pt,fc,cc,kc,CameraModel)

if strcmp(CameraModel,'Tango')
   H_pixel = StaticPixelToHomogeneousDerivativeTango(undistorted_pt,fc,kc);
elseif strcmp(CameraModel,'Radial')
   H_pixel = StaticPixelToHomogeneousDerivativeRadial(undistorted_pt,fc,kc);
elseif strcmp(CameraModel,'Fisheye')
   H_pixel = StaticPixelToHomogeneousDerivativeFisheye(undistorted_pt,fc,kc);
end