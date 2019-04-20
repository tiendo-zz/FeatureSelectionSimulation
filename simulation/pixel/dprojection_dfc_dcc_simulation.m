function [Jfc,Jcc,residual,flag] = dprojection_dfc_dcc_simulation(Cr_p_f, Ci_R_G, G_R_Cr,G_p_Cr, G_p_Ci, Ci_z_f,CameraParams,thresh)

% dfc | dcc

featureNums = size(Cr_p_f,2); 

Ci_p_f = zeros(3,featureNums);
Ci_z_f_hat = zeros(3,featureNums);
Jfc = zeros(2*featureNums, 2);
Jcc = zeros(2*featureNums, 2);
% residual = reshape(Ci_z_f, [2*featureNums 1]) - reshape(Ci_z_f_hat, [2*featureNums 1]);
residual = zeros(2*featureNums, 1);
row_count = 0;
flag = 0;
for k = 0:featureNums-1 
    % transform Cr_p_f to Cartisan
    Cr_p_f_xyz = [cos(Cr_p_f(2,k+1))*cos(Cr_p_f(1,k+1));cos(Cr_p_f(2,k+1))*sin(Cr_p_f(1,k+1));sin(Cr_p_f(2,k+1))];
    Ci_p_f(:,k+1) = Ci_R_G*G_R_Cr*Cr_p_f_xyz + Cr_p_f(3,k+1)*Ci_R_G*(G_p_Cr-G_p_Ci);
    
    homo = Ci_p_f(1:2,k+1) / Ci_p_f(3,k+1);
    Ci_z_f_hat(:,k+1) = CameraParams.K*[homo;1];
    
    % residual
    residual_tmp = Ci_z_f(1:2,k+1) - Ci_z_f_hat(1:2,k+1);
    gamma = norm(residual_tmp);
    if  gamma > thresh % 16.27, chi-square for d.o.f = 3, p = 0.001
        if gamma > thresh*1.5
            flag = 2;
        else
            flag = 1;
        end
        continue;
    end
    
    % Jacobian
    residual(2*row_count+1:2*row_count+2,1) = residual_tmp;
    Jfc(2*row_count+1:2*row_count+2,:) = [homo(1) 0;0 homo(2)];  
    Jcc(2*row_count+1:2*row_count+2,:) = eye(2);
    row_count = row_count + 1;
end

% clear the zero lines

Jfc(2*row_count+1:end,:) = [];
Jcc(2*row_count+1:end,:) = [];
residual(2*row_count+1:end,:) = [];

if row_count == 0
    Jfc = [];
    Jcc = [];
    residual = [];
end