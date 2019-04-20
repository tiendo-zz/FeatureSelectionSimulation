function H_pixel = StaticPixelToHomogeneousDerivativeRadial(undistorted_pt,fc,kc) 
  H_pixel = zeros(2,2);
  

  u = undistorted_pt(1);
  v = undistorted_pt(2);
  u_sq = u * u;
  v_sq = v * v;
  r_sq = u_sq + v_sq;
  r_four = r_sq * r_sq;
  pix_term1 = 2.0 * kc(1) + 4.0 * kc(2) * r_sq ...
      + 6.0 * kc(5) * r_four;
  pix_term2 = 1.0 + kc(1) * r_sq + kc(2) * r_four ...
        + kc(5) * r_four * r_sq + 2.0 * kc(3) * v + 2.0 * kc(4) * u;
  H_pixel(1, 1) = pix_term1 * u_sq + pix_term2 + 4.0 * kc(4) * u;
  H_pixel(2, 2) = pix_term1 * v_sq + pix_term2 + 4.0 * kc(3) * v;
  H_pixel(1, 2) = pix_term1 * u * v + 2.0 * kc(3) * u + 2.0 * kc(4) * v;
  H_pixel(2, 1) = H_pixel(1, 2);
  H_pixel(1,:) = H_pixel(1,:) * fc(1);
  H_pixel(2,:) = H_pixel(2,:) * fc(2);