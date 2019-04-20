function J = Img_Jacobian(c2_p_f)

pix = [c2_p_f(1)/c2_p_f(3); c2_p_f(2)/c2_p_f(3)];
J = 1/c2_p_f(3)*[eye(2) -pix];
