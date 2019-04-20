function homo = StaticUndistortFisheye(fc,cc,kc,z)

N = size(z,2);
homo = zeros(2,N);
for k = 1:N
    
kNumIterations = 20;

u = (z(1,k) - cc(1)) / fc(1);
v = (z(2,k) - cc(2)) / fc(2);
k1 = kc(1);
k2 = kc(2);
k3 = kc(3);
k4 = kc(4);
theta_d = sqrt(u * u + v * v);
theta = theta_d;
% if (std::abs(theta_d) > kAngleTolerance) 
%     // This is the tolerance before theta goes to +/- PI/2
%    return false;
 
for i = 0 : kNumIterations-1
    theta2 = theta * theta;
    theta4 = theta2 * theta2;
    theta6 = theta2 * theta4;
    theta8 = theta4 * theta4;
    gx = theta_d - theta * (1.0 + k1 * theta2 + k2 * theta4 + k3 * theta6 + k4 * theta8);
    dgx = -(1.0 + 3 * k1 * theta2 + 5 * k2 * theta4 + 7 * k3 * theta6 + 9 * k4 * theta8);
    theta = theta - gx / dgx;
end
scaling = tan(theta) / theta_d;
homo(1,k) = u * scaling;
homo(2,k) = v * scaling;

end
