function [AbsolutePoses,PairWiseMatches_homo,c1_f_hat] = Two_View_BLS(PairWiseMatches_homo,AbsolutePoses)

data = struct;
data.C2_q_C1_init = rot2quat(AbsolutePoses(:,1:3));
data.C2_t_C1_init = AbsolutePoses(:,4);
N = size(PairWiseMatches_homo,1);
data.C2_homo = [PairWiseMatches_homo(:,4:5)';ones(1,N)];
data.C1_homo = [PairWiseMatches_homo(:,1:2)';ones(1,N)];

max_iter = 10;
[c2_q_c1, c2_t_c1, c1_f_hat,inlier_id] = TwoViewsReconstruction(data, max_iter);

AbsolutePoses(:,1:3) = quat2rot(c2_q_c1);
AbsolutePoses(:,4) = c2_t_c1;
PairWiseMatches_homo = PairWiseMatches_homo(inlier_id,:);
c1_f_hat = reshape(c1_f_hat,[3,size(c1_f_hat,1)/3]);