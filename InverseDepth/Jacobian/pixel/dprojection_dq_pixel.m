function [Jq,residual,flag] = dprojection_dq(Cr_p_f, Ci_R_G, G_R_Cr,G_p_Cr, G_p_Ci, Ci_z_f)

% dCi_q_G 

featureNums = size(Cr_p_f,2); 


Ci_p_f = zeros(3,featureNums);
Ci_z_f_hat = zeros(2,featureNums);
Jq = zeros(2*featureNums, 3);

residual = zeros(2*featureNums, 1);
row_count = 0;
flag = 0;
for k = 0:featureNums-1
    % if residual of current measurement greater than threshold, erase this
    % measurement.
    
    % transform Cr_p_f to Cartisan
    m_theta_phi = [cos(Cr_p_f(2,k+1))*cos(Cr_p_f(1,k+1));cos(Cr_p_f(2,k+1))*sin(Cr_p_f(1,k+1));sin(Cr_p_f(2,k+1))];
    Ci_p_f(:,k+1) = Ci_R_G*G_R_Cr*m_theta_phi + Cr_p_f(3,k+1)*Ci_R_G*(G_p_Cr-G_p_Ci);
    Ci_z_f_hat(:,k+1) = Ci_p_f(1:2,k+1) / Ci_p_f(3,k+1);
    
    residual_tmp = Ci_z_f(1:2,k+1) - Ci_z_f_hat(1:2,k+1);
    gamma = residual_tmp'*10^4*residual_tmp;
    if  gamma > 16.27 % 16.27, chi-square for d.o.f = 3, p = 0.001
        if gamma > 2*16.27
            flag = 2;
        else
            flag = 1;
        end
        continue;
    end
    residual(2*row_count+1:2*row_count+2,1) = residual_tmp;
    Jq(2*row_count+1:2*row_count+2,1:3) = 1/Ci_p_f(3,k+1)*[1 0 -Ci_p_f(1,k+1)/Ci_p_f(3,k+1); ...
                                            0 1 -Ci_p_f(2,k+1)/Ci_p_f(3,k+1)]*skewsymm(Ci_p_f(:,k+1));
    row_count = row_count + 1;
end

% clear the zero lines
Jq(2*row_count+1:end,:) = [];
residual(2*row_count+1:end,:) = [];

if row_count == 0
    Jq = [];
    residual = [];
end