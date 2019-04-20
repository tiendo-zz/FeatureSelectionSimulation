function H_pixel = StaticPixelToHomogeneousDerivativeTango(undistorted_pt,fc,kc)

  u = undistorted_pt(1);
  v = undistorted_pt(2);
  h_x_sq = u * u;
  h_y_sq = v * v;
  r_sq = h_x_sq + h_y_sq;
  r = sqrt(r_sq);
  r_cub = r_sq * r;
  r_four = r_sq * r_sq;
  tankc = tan(kc(1) / 2);
  beta = atan(2 * r * tankc);
  H_pixel(1, 1) = fc(1) / kc(1) ...
        * (beta / r + 2 * h_x_sq * tankc / (r_sq + 4 * r_four * tankc * tankc)...
            - h_x_sq * beta / r_cub);
  H_pixel(2, 2) = fc(2) / kc(1) ...
        * (beta / r + 2 * h_y_sq * tankc / (r_sq + 4 * r_four * tankc * tankc)...
            - h_y_sq * beta / r_cub);
  cross_term = 1.0 / kc(1) * u * v ...
        * (2 * tankc / (r_sq + 4 * r_four * tankc * tankc) - beta / r_cub);
  H_pixel(2, 1) = fc(2) * cross_term;
  H_pixel(1, 2) = fc(1) * cross_term;
end

