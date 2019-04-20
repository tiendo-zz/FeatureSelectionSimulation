function  pixel = StaticDistortFisheye(fc,cc,kc,point)

N = size(point,2);
pixel = zeros(2,N);


if size(cc,2) > 1
    cc= cc';
end
if size(fc,2) > 1
    fc= fc';
end

for i = 1:N
	% Convert from pixels to camera center frame
	point_normalized = point(:,i) - cc;
	point_normalized(1) = point_normalized(1) / fc(1);
	point_normalized(2) = point_normalized(2) / fc(2);

	r = norm(point_normalized);

	theta = atan(r);

	k1 = kc(1);
	k2 = kc(2);
	k3 = kc(3);
	k4 = kc(4);

	theta2 = theta * theta;
	theta4 = theta2 * theta2;
	theta6 = theta2 * theta4;
	theta8 = theta4 * theta4;

	theta_d = theta * (1 + k1 * theta2 + k2 * theta4 + k3 * theta6 + k4 * theta8);

    if (r > 1e-8)
        scaling = theta_d / r;
    else
        scaling = 1;
    end
    
	% compute distorted in camera center frame
	distorted = scaling*point_normalized;

	% convert from camera center to pixels
	distorted(1) = distorted(1)*fc(1);
	distorted(2) = distorted(2)*fc(2);
	pixel(:,i) = distorted + cc;
end
