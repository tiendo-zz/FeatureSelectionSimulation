function  z = StaticDistortTango(fc,cc,kc,homo)

% change the homo value to pixel value
N = size(homo,2);
z = zeros(2,N);

for i = 1:N
u = homo(1,i);
v = homo(2,i);
r = (u^2 + v^2)^0.5;
distortion_factor = atan(2 * r * tan(kc(1) / 2)) / (kc(1) * r);
x_d = fc(1) * distortion_factor * u + cc(1);
y_d = fc(2) * distortion_factor * v + cc(2);
z(:,i) = [x_d;y_d];
end
