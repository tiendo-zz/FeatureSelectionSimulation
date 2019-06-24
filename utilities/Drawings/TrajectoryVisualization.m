% Visualization script for trajectory
NumOfPoses = 10;Rcam=1;
angle_range = 2*pi; WpC_ = [];
figure(1); hold on;
for k = 1:NumOfPoses
    tt = angle_range*(k-1)/NumOfPoses;
    phi = angle_range*(k-1)/NumOfPoses;  
    WpC = Rcam*[cos(tt)*cos(phi);sin(tt)*cos(phi);sin(phi)];
    WpC_ = [WpC_;WpC'];
    k_axis = skewsymm(WpC)*[0;0;-1]; 
    k_axis = k_axis/norm(k_axis);
    CRW = aa2rot(k_axis, acos(-WpC(3)/norm(WpC)));
    WRC = CRW'
    
    
    s = sin(tt); c = cos(tt);
    R = [-s+(1+s)*s^2 (1+s)*(-s*c) -c^2;...
         (1+s)*(-s*c) -s+(1+s)*c^2 -s*c;...
         c^2 c*s -s]
    
    
    if k >= 2
    s = plot3(WpC_(end-1:end,1), WpC_(end-1:end,2), WpC_(end-1:end,3), 'b-'); grid on;    
    set(s, 'LineWidth',3)
    axis([-1 1 -1 1 -1 1]);
    pause(0.02);
    end
end