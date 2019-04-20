function residual_r = Robust_CostFunction(residual,sigma)

if norm(residual)^2 < 1e-6
    x_r = 1/sigma*residual(1);
    y_r = 1/sigma*residual(2);
else
    x_r = (log(1+norm(residual)^2/sigma^2)/norm(residual)^2)^0.5*residual(1);
    y_r = (log(1+norm(residual)^2/sigma^2)/norm(residual)^2)^0.5*residual(2);
end
residual_r = sigma*[x_r;y_r];