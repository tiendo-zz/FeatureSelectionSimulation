function Jr = Robust_Cost_Jacobian(residual,sigma)

x = residual(1);
y = residual(2);

if norm(residual)^2 < 1e-4
    f_xy = 1/sigma;
    if y/x > 100
        y_devide_norm_xy =  1;
        x_devide_norm_xy =  0;
    elseif x/y > 100
        y_devide_norm_xy =  0;
        x_devide_norm_xy =  1;
    else
        y_devide_norm_xy =  y^2/(x^2+y^2);
        x_devide_norm_xy =  x^2/(x^2+y^2);
    end
        
else
    f_xy = (log(1+(x^2+y^2)/sigma^2)/(x^2+y^2))^0.5;
    y_devide_norm_xy =  y^2/(x^2+y^2);
    x_devide_norm_xy =  x^2/(x^2+y^2);
end


Jxx = y_devide_norm_xy*f_xy + x_devide_norm_xy*1/(sigma^2+x^2+y^2)*1/f_xy;
Jxy = x*y/(x^2+y^2)*(1/(f_xy*(sigma^2+x^2+y^2))-f_xy);
Jyy = x_devide_norm_xy*f_xy + y_devide_norm_xy*1/(sigma^2+x^2+y^2)*1/f_xy;

Jr = sigma*[Jxx Jxy;Jxy Jyy];