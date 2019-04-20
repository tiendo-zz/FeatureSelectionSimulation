function homo = StaticUndistortRadial(fc,cc,kc,pixel)

N = size(pixel,2);
homo = zeros(2,N);

for i = 1:N    
   homo(1,i) = (pixel(1,i) - cc(1)) / fc(1);
   homo(2,i) = (pixel(2,i) - cc(2)) / fc(2);
 
   k1 = kc(1);
   k2 = kc(2);
   k5 = kc(5);
   k3 = kc(3);
   k4 = kc(4);
 
   k_radial = 0;
   r_2 = 0;
   xd0 = homo(1,i);
   xd1 = homo(2,i);
   
   for k = 1:20
     p00 = homo(1,i) * homo(1,i);  %point(0) * point(0);
     p11 = homo(2,i) * homo(2,i);  %point(1) * point(1);
     p01 = homo(1,i) * homo(2,i);  %point(0) * point(1);
     r_2 = p00 + p11;
     k_radial = 1 / (1 + (k1 * r_2) + (k2 * r_2 * r_2) + (k5 * r_2 * r_2 * r_2));
     

     A=[6*k3*homo(1,i) + 2*k4*homo(2,i), 2*k3*homo(2,i) + 2*k4*homo(1,i);
        2*k4*homo(2,i) + 2*k3*homo(1,i), 2*k3*homo(2,i) + 2*k4*homo(1,i)];
 
     B = eye(2);
     B = B - k_radial*A;
 
     delta = zeros(2,1);
     delta(1) = (xd0 - (2*k3*p01 + k4*(r_2+2*p00)));
     delta(2) = (xd1 - (k3*(r_2+2*p11) + 2*k4*p01));
 
     homo(:,i) = A*homo(:,i) + B*delta;
     homo(:,i) = homo(:,i)*k_radial;
   end
    
end

end