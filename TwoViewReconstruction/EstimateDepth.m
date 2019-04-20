function depths = EstimateDepth(c2_b_f, c1_b_f, c2_q_c1, c2_t_c1)

A = [-quat2rot(c2_q_c1)*c1_b_f c2_b_f];
depths = (A'*A) \ (A'*c2_t_c1);
