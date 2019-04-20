function [v_p, v_pp] = Orthonormal_3D_Set(v)

v = v/norm(v);
[~, i] = max(abs(v));
c = v(i);
v(i) = 0;
[~, m] = max(abs(v));
v_p = zeros(3,1);
v_p(m) = c;
v_p(i) = -v(m);
v(i) = c;
v_pp = cross(v, v_p);

v_p = v_p/norm(v_p);
v_pp = v_pp/norm(v_pp);
