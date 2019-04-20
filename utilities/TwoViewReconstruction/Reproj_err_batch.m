function [err, c2_f_hat] = Reproj_err_batch(c2_q_c1, c2_t_c1, c1_p_f, c2_homo, c1_homo)

N = size(c2_homo,2);
c2_R_c1 = quat2rot(c2_q_c1);
err = 0;
c2_f_hat = [];
for i = 1:N
    c1_p_f_i = c1_p_f(3*(i-1)+1:3*i);
    c2_p_f_i = c2_R_c1 * c1_p_f_i + c2_t_c1;
    c2_f_hat = [c2_f_hat; c2_p_f_i];
    % error in first view
    err = err + Reproj_err(c1_homo(:,i), c1_p_f_i);
    
    % error in second view
    err = err + Reproj_err(c2_homo(:,i), c2_p_f_i);
end

