clear all;
close all;
clc

residual_true = -1.5e-2 + 3e-2*rand(100,1);

sigma_r = 1e-3;
pertube = sigma_r.*rand(100,1);
residual = residual_true - pertube;

error = zeros(length(residual)/2,2);
measure_count = 1;
sigma = 1e-2;
for i = 1:length(residual)/2
    residual_single = residual(2*i-1:2*i);
    residual_r = Robust_CostFunction(residual_single,sigma);
    Jr = Robust_Cost_Jacobian(residual_single,sigma)
    res_ = residual_r + Jr*pertube(2*i-1:2*i);
           
    res_r_true = Robust_CostFunction(residual_true(2*i-1:2*i),sigma);
    % error
    error(measure_count,1) = norm(res_ - res_r_true);
    error(measure_count,2) = norm(residual_r - res_r_true);
    measure_count = measure_count + 1;
end

figure;
percent = sum(error(:,1)<error(:,2))/size(error,1)
plot(error(:,1));hold on
plot(error(:,2));
legend('measurement ','estimation');