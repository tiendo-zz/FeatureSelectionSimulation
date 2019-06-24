function [time, trajectory, measurements] = real_world_simulation(time_config, trajectory_config, measurements_config)
% Real world simulation: Generate data for a slam experiment, according to
% the specified configuration.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate timing data

% Count of experiment timesteps
time.count_time_steps = round(time_config.duration / time_config.time_step) + 1;

% True time steps
time.time_step_true = [0, time_config.time_step + time_config.sigma_time_step * randn(1, time.count_time_steps - 1)];

% True absolute time
time.absolute_time_true = cumsum(time.time_step_true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate true trajectory data

% True orientation
trajectory.position_true = zeros(3, time.count_time_steps);
for i = 1:time.count_time_steps
    trajectory.position_true(:,i) = trajectory_config.position(time.absolute_time_true(i));
end

% True orientation
trajectory.orientation_true = zeros(3, 3, time.count_time_steps);
for i = 1:time.count_time_steps
    trajectory.orientation_true(:,:,i) = trajectory_config.orientation(time.absolute_time_true(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate measurement data

% True relative 6-DoF pose (orientation, position).
measurements.relative_pose_true = zeros(7, time.count_time_steps);

for i = 2 : time.count_time_steps
    measurements.relative_pose_true(1:4, i) = rotm2quat(trajectory.orientation_true(:,:,i-1)'*trajectory.orientation_true(:,:,i))';
    measurements.relative_pose_true(5:7, i) = trajectory.orientation_true(:,:,i-1)'*(trajectory.position_true(:,i) - trajectory.position_true(:,i-1));
end

% Relative 6-DoF pose (orientation, position) measurements.
error_relative_pose = chol(measurements_config.cov_relative_pose,'lower') * randn(6,time.count_time_steps); 
error_relative_orientation = [sqrt(1 - 0.25*vecnorm(error_relative_pose(1:3,:)).^2)', 0.5*error_relative_pose(1:3,:)'];
error_relative_position = error_relative_pose(4:6,:);
measurements.relative_pose_meas = [quatmultiply(error_relative_orientation, measurements.relative_pose_true(1:4, :)')'; 
                                                     measurements.relative_pose_true(5:7, :) + error_relative_position];

end
