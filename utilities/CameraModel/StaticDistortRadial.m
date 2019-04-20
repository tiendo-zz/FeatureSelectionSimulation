function pixel = StaticDistortRadial(fc,cc,kc,homo)

N = size(homo,2);
pixel = zeros(2,N);

for i = 1:N
  u = homo(1,i);
  v = homo(2,i);
  r_sq = u * u + v * v;
  r_four = r_sq * r_sq;
  pix_term2 = 1.0 + kc(1) * r_sq + kc(2) * r_four + kc(5) * r_four * r_sq + 2.0 * kc(3) * v + 2.0 * kc(4) * u;
  pixel(1,i) = fc(1) * (pix_term2 * u + kc(4) * r_sq) + cc(1);
  pixel(2,i) = fc(2) * (pix_term2 * v + kc(3) * r_sq) + cc(2);
end