function homo = StaticUndistortTango(fc,cc,kc,z)

N = size(z,2);
homo = zeros(2,N);

for i = 1:N
    kAngleLimit = pi * (89.0 / 180.0);
    u = (z(1,i) - cc(1)) / fc(1);
    v = (z(2,i) - cc(2)) / fc(2);
    omega = kc(1); 
    r_d = sqrt(u * u + v * v);
    mul2tanwby2 = tan(omega / 2) * 2;

    if (abs(r_d * mul2tanwby2) < 1e-6)
         homo(:,i) = [u;v];
    end

    if (abs(r_d * omega) <= kAngleLimit) 
        r_u = tan(r_d * omega) / (r_d * mul2tanwby2);
        homo(1,i) = u * r_u;
        homo(2,i) = v * r_u;
    end
end
