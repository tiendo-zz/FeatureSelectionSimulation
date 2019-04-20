function Jtheta = dprojection_dtheta(Cr_p_f, Ci_R_G, G_R_Cr,G_p_Cr, G_p_Ci, Ci_z_f)

% dG_theta_Cr

featureNums = size(Cr_p_f,2); 

Ci_p_f = zeros(3,featureNums);
Ci_z_f_hat = zeros(2,featureNums);
Jtheta = zeros(2*featureNums, 1);

% residual = reshape(Ci_z_f, [2*featureNums 1]) - reshape(Ci_z_f_hat, [2*featureNums 1]);
residual = zeros(2*featureNums, 1);
row_count = 0;
for k = 0:featureNums-1
    % if residual of current measurement greater than threshold, erase this
    % measurement.
    
     % transform Cr_p_f to Cartisan
    Cr_p_f_xyz = [cos(Cr_p_f(2,k+1))*cos(Cr_p_f(1,k+1));cos(Cr_p_f(2,k+1))*sin(Cr_p_f(1,k+1));sin(Cr_p_f(2,k+1))];
    Ci_p_f(:,k+1) = Ci_R_G*G_R_Cr*Cr_p_f_xyz + Cr_p_f(3,k+1)*Ci_R_G*(G_p_Cr-G_p_Ci);
    Ci_z_f_hat(:,k+1) = Ci_p_f(1:2,k+1) / Ci_p_f(3,k+1);
    
    residual_tmp = Ci_z_f(1:2,k+1) - Ci_z_f_hat(1:2,k+1);
    residual(2*row_count+1:2*row_count+2,1) = residual_tmp;
    dm_dtheta = [-cos(Cr_p_f(2,k+1))*sin(Cr_p_f(1,k+1));cos(Cr_p_f(2,k+1))*cos(Cr_p_f(1,k+1));0];
    Jtheta(2*row_count+1:2*row_count+2,:) = 1/Ci_p_f(3,k+1)*[1 0 -Ci_p_f(1,k+1)/Ci_p_f(3,k+1); ...
                                            0 1 -Ci_p_f(2,k+1)/Ci_p_f(3,k+1)]*(Ci_R_G*G_R_Cr*dm_dtheta);
    row_count = row_count + 1;
end
