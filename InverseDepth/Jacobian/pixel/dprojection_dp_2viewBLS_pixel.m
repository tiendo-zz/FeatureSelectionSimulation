function [Jp,residual,flag,t_p, t_pp] = dprojection_dp_2viewBLS_pixel(Cr_p_f, Ci_R_G, G_R_Cr,G_p_Cr, G_p_Ci, Ci_z_f)

% dG_p_Ci

featureNums = size(Cr_p_f,2); 

Ci_p_f = zeros(3,featureNums);
Ci_z_f_hat = zeros(2,featureNums);
Jp = zeros(2*featureNums, 2);

residual = zeros(2*featureNums, 1);
row_count = 0;
flag = 0;
[t_p, t_pp] = Orthonormal_3D_Set(G_p_Ci);
for k = 0:featureNums-1
    % if residual of current measurement greater than threshold, erase this
    % measurement.
    
    % transform Cr_p_f to Cartisan
    Cr_p_f_xyz = [cos(Cr_p_f(2,k+1))*cos(Cr_p_f(1,k+1));cos(Cr_p_f(2,k+1))*sin(Cr_p_f(1,k+1));sin(Cr_p_f(2,k+1))];
    Ci_p_f(:,k+1) = Ci_R_G*G_R_Cr*Cr_p_f_xyz + Cr_p_f(3,k+1)*Ci_R_G*(G_p_Cr-G_p_Ci);
    Ci_z_f_hat(:,k+1) = Ci_p_f(1:2,k+1) / Ci_p_f(3,k+1);
    
    residual_tmp = Ci_z_f(1:2,k+1) - Ci_z_f_hat(1:2,k+1);
%     gamma = residual_tmp'*10^4*residual_tmp;
%     if  gamma > 16.27  % chi-square for d.o.f = 3, p = 0.001
%         if gamma > 2*16.27
%             flag = 2;
%         else
%             flag = 1;
%         end
%         continue;
%     end
   residual(2*row_count+1:2*row_count+2,1) = residual_tmp;
   
   Jp(2*row_count+1:2*row_count+2,:) = 1/Ci_p_f(3,k+1)*[1 0 -Ci_p_f(1,k+1)/Ci_p_f(3,k+1); ...
                                            0 1 -Ci_p_f(2,k+1)/Ci_p_f(3,k+1)]*(-Cr_p_f(3,k+1)*Ci_R_G)*[-t_p t_pp];
   row_count = row_count + 1;
end

% clear the zero lines
Jp(2*row_count+1:end,:) = [];
residual(2*row_count+1:end,:) = [];

if row_count == 0
    Jp = [];
    residual = [];
end