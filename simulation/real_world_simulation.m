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
error_relative_orientation = [sqrt(1 - 0.25*vecnorm(dqp(1:3,:)).^2)', 0.5*dqp(1:3,:)'];
error_relative_position = error_relative_pose(4:6,:);



% Relative 6DoF pose measurements availability
ricp_meas_avail = false(1,N);
ricp_meas_avail(ricp_rate+1:ricp_rate:N) = true;
ricp_meas_idx_info = zeros(3,N);
ricp_meas_idx_info(1,ricp_rate+1:ricp_rate:N) = 1;
ricp_meas_idx_info(2,ricp_rate+1:ricp_rate:N) = 1:ricp_rate:N-ricp_rate;
ricp_meas_idx_info(3,ricp_rate+1:ricp_rate:N) = ricp_rate+1:ricp_rate:N;

ricp_true = zeros(7,N);
for i = find(ricp_meas_avail)
    idx_meas = i;
    idx_start = ricp_meas_idx_info(2,i);
    idx_end = ricp_meas_idx_info(3,i);
    ricp_true(1:4,idx_meas) = rotm2quat(G_R_B(t(idx_start))'* G_R_B(t(idx_end)));
    ricp_true(5:7,idx_meas) = G_R_B(t(idx_start))'*(G_p_B(t(idx_end)) - G_p_B(t(idx_start)));   
end

idx_meas = find(ricp_meas_idx_info(1,:));
dqp = chol(cov_ricp,'lower') * randn(6,length(idx_meas)); 
dq = [sqrt(1 - 0.25*vecnorm(dqp(1:3,:)).^2)', 0.5*dqp(1:3,:)'];
dp = dqp(4:6,:);

ricp_meas = zeros(7,N);
ricp_meas(:,idx_meas) = [quatmultiply(dq, ricp_true(1:4,idx_meas)')'; ricp_true(5:7,idx_meas) + dp];






end
