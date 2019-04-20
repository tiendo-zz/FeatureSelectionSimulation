function fc_jacobian = StaticPixelToCameraIntrinsicsDerivative(undistorted_pt,kc,CameraModel)

fc_jacobian = zeros(2,2);
if strcmp(CameraModel,'Tango')
    u = undistorted_pt(1);
    v = undistorted_pt(2);
    r = sqrt(u * u + v * v);
    distortion_factor = atan(2 * r * tan(kc(1) / 2)) / (kc(1) * r);
    fc_jacobian(1, 1) = distortion_factor * u;
    fc_jacobian(1, 2) = 0;
    fc_jacobian(2, 1) = 0;
    fc_jacobian(2, 2) = distortion_factor * v;
elseif strcmp(CameraModel,'Radial')
    u = undistorted_pt(1);
    v = undistorted_pt(2);
    u_sq = u * u;
    v_sq = v * v;
    r_sq = u_sq + v_sq;
    r_four = r_sq * r_sq;
    pix_term2 = 1.0 + kc(1) * r_sq + kc(2) * r_four ...
          + kc(5) * r_four * r_sq + 2.0 * kc(3) * v + 2.0 * kc(4) * u;
    fc_jacobian(1, 1) = pix_term2 * u + kc(4) * r_sq;
    fc_jacobian(1, 2) = 0;
    fc_jacobian(2, 1) = 0;
    fc_jacobian(2, 2) = pix_term2 * v + kc(3) * r_sq;
end