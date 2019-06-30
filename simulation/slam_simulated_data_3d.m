%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load dependencies
addpath('../utilities/GeometricToolbox/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load configuration settings
% TODO: Move implementation of trajectory generators in a function

% Timing configuration parameters
time_config.duration = 20;                 % Duration of the experiment (s)
time_config.time_step = 1;                 % System timestep (s)
time_config.sigma_time_step = 0.0001;      % Timestep gitter profile (s)

% Trajectory configuration parameters
trajectory_config.group = 'infinity_morbius';
trajectory_config.radius = 35;            % Trajectory sphere radius (m)
trajectory_config.angular_freq = 2*pi/20;  % Trajectory angular frequency (rad/s)
trajectory_config.position = @(t)...
    trajectory_config.radius * [cos(trajectory_config.angular_freq*t)*cos(trajectory_config.angular_freq*t);
                                sin(trajectory_config.angular_freq*t)*cos(trajectory_config.angular_freq*t);
                                sin(trajectory_config.angular_freq*t);];
trajectory_config.orientation = @(t)...
    [-sin(trajectory_config.angular_freq*t)+(1+sin(trajectory_config.angular_freq*t))*sin(trajectory_config.angular_freq*t)^2, (1+sin(trajectory_config.angular_freq*t))*(-sin(trajectory_config.angular_freq*t)*cos(trajectory_config.angular_freq*t)),                                     -cos(trajectory_config.angular_freq*t)^2;...
     (1+sin(trajectory_config.angular_freq*t))*(-sin(trajectory_config.angular_freq*t)*cos(trajectory_config.angular_freq*t)), -sin(trajectory_config.angular_freq*t)+(1+sin(trajectory_config.angular_freq*t))*cos(trajectory_config.angular_freq*t)^2, -sin(trajectory_config.angular_freq*t)*cos(trajectory_config.angular_freq*t);...
                                                                                      cos(trajectory_config.angular_freq*t)^2,                                              cos(trajectory_config.angular_freq*t)*sin(trajectory_config.angular_freq*t),                                       -sin(trajectory_config.angular_freq*t)];
% Measurements configuration parameters
% 6-dof relative pose (orientation, position) measurement parameters 
% measurements_config.cov_relative_pose = [0.01^2, 0, 0,     0, 0, 0;
%                                          0, 0.01^2, 0,     0, 0, 0;
%                                          0, 0, 0.01^2,     0, 0, 0;
%                                          0, 0, 0,      0.1^2, 0, 0;
%                                          0, 0, 0,      0, 0.1^2, 0;
%                                          0, 0, 0,      0, 0, 0.1^2];   
measurements_config.cov_relative_pose = 1e-3 * eye(6);
                                     
delay_line = 20;
                             
% World configuration paprameters
camera_config.fov = pi/3;
camera_config.fc = 400;
camera_config.cc = [640; 480]/2;
camera_config.sigma_pixel = 10; %(pixel)


world_config.group = 'sphere';
world_config.radius = 30;
world_config.center = [0;0;0];
world_config.point_feature_num = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate and plot/visualize data for a slam experiment

% Real world simulation: Generate data for a slam experiment
[time, trajectory, world, measurements] = real_world_simulation(time_config, ...
                                                                trajectory_config, ...
                                                                world_config, ...
                                                                camera_config, ...
                                                                measurements_config);

% Visualize the trajectory and enviroment
% TODO: Make the trajectory visualization into animation
% TODO: Make the ploting for 
figure('Name','True trajectory'); hold on;
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
grid on
axis equal 
for i=1:time.count_time_steps
    p_true = trajectory.position_true(:,i);
    R_true = trajectory.orientation_true(:,:,i);
    scatter3(p_true(1), p_true(2), p_true(3), 15, 'g' , 'filled');
    quiver3(p_true(1), p_true(2), p_true(3),...
            R_true(1,1), R_true(2,1), R_true(3,1),...
            'color', 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver3(p_true(1), p_true(2), p_true(3),...
            R_true(1,2), R_true(2,2), R_true(3,2),...
            'color', 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver3(p_true(1), p_true(2), p_true(3),...
            R_true(1,3), R_true(2,3), R_true(3,3),...
            'color', 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
end
for k = 1:world_config.point_feature_num
    scatter3(world.point_feat_Wposition{k}(1), ...
             world.point_feat_Wposition{k}(2), ...
             world.point_feat_Wposition{k}(3), ...
             'b', 'filled');
end

% Plot the true and measured relative 6DoF poses
relative_rpy_true = flipud(quat2eul(measurements.relative_pose_true(1:4,:)','ZYX')');
relative_rpy_meas = flipud(quat2eul(measurements.relative_pose_meas(1:4,:)','ZYX')');
alpha = 2*pi*round((relative_rpy_true - relative_rpy_meas)/(2*pi));
relative_rpy_meas = relative_rpy_meas + alpha;
error_relative_rpy = zeros(3, time.count_time_steps);
for i=1:size(error_relative_rpy,2)
    roll = relative_rpy_true(1,i);
    pitch = relative_rpy_true(2,i);
    yaw = relative_rpy_true(3,i);
    H_rpy = [cos(yaw)/cos(pitch), sin(yaw)/cos(pitch), 0;
                       -sin(yaw),            cos(yaw), 0;
             cos(yaw)*tan(pitch), sin(yaw)*tan(pitch), 1];        
    error_relative_rpy(:,i) = sqrt(diag(H_rpy * measurements_config.cov_relative_pose(1:3,1:3) * H_rpy')); 
end
figure('Name','Relative 6DoF orientation'); hold on
angle = {'roll (rad)','pitch (rad)','yaw (rad)'};
for i = 1:3
    subplot(3,1,i);
    plot(time.absolute_time_true, relative_rpy_meas(i,:), 'b'); hold on;
    plot(time.absolute_time_true, relative_rpy_true(i,:), 'g'); hold on;
    plot(time.absolute_time_true, relative_rpy_true(i,:) + 3*error_relative_rpy(i,:),'r'); hold on;
    plot(time.absolute_time_true, relative_rpy_true(i,:) - 3*error_relative_rpy(i,:),'r'); hold on;
    xlabel('time (s)')
    ylabel(angle{i})
    legend('measured','true','3\sigma - bound')
    grid on
end
error_relative_position = sqrt(diag(measurements_config.cov_relative_pose(4:6,4:6)));
figure('Name','Relative 6DoF position'); hold on
axes = {'x (m)','y (m)','z (m)'};
for i = 1:3
    subplot(3,1,i);
    plot(time.absolute_time_true, measurements.relative_pose_meas(4+i,:), 'b'); hold on;
    plot(time.absolute_time_true, measurements.relative_pose_true(4+i,:),'g'); hold on;
    plot(time.absolute_time_true, measurements.relative_pose_true(4+i,:) + 3*error_relative_position(i),'r'); hold on;
    plot(time.absolute_time_true, measurements.relative_pose_true(4+i,:) - 3*error_relative_position(i),'r'); hold on;
    xlabel('time (s)')
    ylabel(axes{i})
    legend('measured','true','3\sigma - bound')
    grid on
end        

% scattering corner measurements of a random camera
figure('Name','Example corner measurements'); hold on;
for k = 5
    scatter(measurements.corner_pixel_true{k}(1,:), ...
            measurements.corner_pixel_true{k}(2,:), ...
             'bo');
    scatter(measurements.corner_pixel{k}(1,:), ...
            measurements.corner_pixel{k}(2,:), ...
            'rx');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate and plot/visualize data for a slam experiment

ego_state.size = 7;
ego_state.ori_size = 4; 
ego_state.ori_idx = 1:4;
ego_state.pos_size = 3; 
ego_state.pos_idx = 5:7;
ego_state.pose_size = ego_state.ori_size + ego_state.pos_size;
ego_state.pose_idx = [ego_state.ori_idx, ego_state.pos_idx];
ego_state.count = 1;

world_state.point_feat_size = 3;
world_state.point_feat_idx = 1:3;
world_state.point_feat_count = world_config.point_feature_num;


error_ego_state.size = 6;
error_ego_state.ori_size = 3; 
error_ego_state.ori_idx = 1:3;
error_ego_state.pos_size = 3; 
error_ego_state.pos_idx = 4:6;
error_ego_state.pose_size = error_ego_state.ori_size + error_ego_state.pos_size;
error_ego_state.pose_idx = [error_ego_state.ori_idx, error_ego_state.pos_idx];
error_ego_state.pose_count = 1;

error_world_state.point_feat_size = 3;
error_world_state.point_feat_idx = 1:3;
error_world_state.point_feat_count = world_config.point_feature_num;

x_est = cell(1, time.count_time_steps);
P_est = cell(1, time.count_time_steps);

% estimated state and covariance
% this script initialize the size and
for i= 1:time.count_time_steps
    window_size = i - max(i-delay_line,1) + 1;
        
    x_est{i}.sz_idx.ego = ego_state;
    x_est{i}.sz_idx.ego.count = window_size;
    x_est{i}.sz_idx.ego_idx = cell(window_size,1);
    % Fill the ego index
    for j=1:x_est{i}.sz_idx.ego.count
        x_est{i}.sz_idx.ego_idx{j} = ...
            ((j-1)*x_est{i}.sz_idx.ego.size+1):(j*x_est{i}.sz_idx.ego.size);
    end
    x_est{i}.sz_idx.world = error_world_state;
    x_est{i}.sz_idx.world_idx = cell(x_est{i}.sz_idx.world.point_feat_count,1);
    % Fill the world index
    for j=1:x_est{i}.sz_idx.world.point_feat_count
        x_est{i}.sz_idx.world_idx{j} = ...
            ((j-1)*x_est{i}.sz_idx.world.point_feat_size+1):(j*x_est{i}.sz_idx.world.point_feat_size);
    end
    x_est{i}.data = zeros(x_est{i}.sz_idx.ego.count*x_est{i}.sz_idx.ego.size + ...
                          x_est{i}.sz_idx.world.point_feat_count*x_est{i}.sz_idx.world.point_feat_size, 1);
    
    
    P_est{i}.sz_idx.ego = error_ego_state;
    P_est{i}.sz_idx.ego.count = window_size;
    P_est{i}.sz_idx.ego_idx = cell(window_size,1);
    % Fill the ego index
    for j=1:P_est{i}.sz_idx.ego.count
        P_est{i}.sz_idx.ego_idx{j}=...
            ((j-1)*P_est{i}.sz_idx.ego.size+1):(j*P_est{i}.sz_idx.ego.size);
    end
    
    P_est{i}.sz_idx.world = world_state;
    P_est{i}.sz_idx.world_idx = cell(x_est{i}.sz_idx.world.point_feat_count,1);    
    % Fill the world index
    for j=1:x_est{i}.sz_idx.world.point_feat_count
        x_est{i}.sz_idx.world_idx{j} = ...
            ((j-1)*x_est{i}.sz_idx.world.point_feat_size+1):(j*x_est{i}.sz_idx.world.point_feat_size);
    end
    full_error_state_size = P_est{i}.sz_idx.ego.count*P_est{i}.sz_idx.ego.size + ...
                            P_est{i}.sz_idx.world.point_feat_count*P_est{i}.sz_idx.world.point_feat_size;
    P_est{i}.data = zeros(full_error_state_size,full_error_state_size);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation only
% TODO: Pack it into a function
for i = 1:time.count_time_steps    
    if i == 1
        % Initialize state
        ori_state_ids = x_est{i}.sz_idx.ego_idx{i}(x_est{i}.sz_idx.ego.ori_idx);
        pos_state_ids = x_est{i}.sz_idx.ego_idx{i}(x_est{i}.sz_idx.ego.pos_idx);
        x_est{i}.data(pos_state_ids) = trajectory.position_true(:,i);
        x_est{i}.data(ori_state_ids) = rotm2quat(trajectory.orientation_true(:,:,i))';
        % Initialize covariance
        error_ori_ids = P_est{i}.sz_idx.ego_idx{i}(P_est{i}.sz_idx.ego.ori_idx);
        error_pos_ids = P_est{i}.sz_idx.ego_idx{i}(P_est{i}.sz_idx.ego.pos_idx);
        P_est{i}.data([error_ori_ids error_pos_ids],[error_ori_ids error_pos_ids]) = ...
            1e-6 * eye(P_est{i}.sz_idx.ego.size);
        continue;
    end
    ori_state_ids = x_est{i}.sz_idx.ego_idx{i}(x_est{i}.sz_idx.ego.ori_idx);
    pos_state_ids = x_est{i}.sz_idx.ego_idx{i}(x_est{i}.sz_idx.ego.pos_idx);
    ori_state_prev_ids = x_est{i-1}.sz_idx.ego_idx{i-1}(x_est{i-1}.sz_idx.ego.ori_idx);
    pos_state_prev_ids = x_est{i-1}.sz_idx.ego_idx{i-1}(x_est{i-1}.sz_idx.ego.pos_idx);
    
    % copy/clone the previous state
    x_est{i}.data(1:x_est{i-1}.sz_idx.ego.count*x_est{i-1}.sz_idx.ego.size) = ...
        x_est{i-1}.data(1:x_est{i-1}.sz_idx.ego.count*x_est{i-1}.sz_idx.ego.size);
    
    % augment new state from rel. pose measurement
    x_est{i}.data(ori_state_ids) = quatmultiply( x_est{i-1}.data(ori_state_prev_ids)', measurements.relative_pose_meas(x_est{i}.sz_idx.ego.ori_idx, i)' )';
    x_est{i}.data(pos_state_ids) = x_est{i-1}.data(pos_state_prev_ids) + ...
        quat2rotm(x_est{i-1}.data(ori_state_prev_ids)')*measurements.relative_pose_meas(x_est{i}.sz_idx.ego.pos_idx, i);        
    
    % cov propagation
    ego_past_state = 1:P_est{i-1}.sz_idx.ego.count*P_est{i-1}.sz_idx.ego.size;
    P_est{i}.data(ego_past_state,ego_past_state) = ...
        P_est{i-1}.data(ego_past_state,ego_past_state);
    error_ori_ids = P_est{i}.sz_idx.ego_idx{i}(P_est{i}.sz_idx.ego.ori_idx);
    error_pos_ids = P_est{i}.sz_idx.ego_idx{i}(P_est{i}.sz_idx.ego.pos_idx);
    error_ori_prev_ids = P_est{i-1}.sz_idx.ego_idx{i-1}(P_est{i-1}.sz_idx.ego.ori_idx);
    error_pos_prev_ids = P_est{i-1}.sz_idx.ego_idx{i-1}(P_est{i-1}.sz_idx.ego.pos_idx);
    Phi = [eye(3) zeros(3);...
        -skewsymm(quat2rotm(x_est{i-1}.data(ori_state_prev_ids)')*measurements.relative_pose_meas(x_est{i}.sz_idx.ego.pos_idx, i)) eye(3)];
    G = -[quat2rotm(x_est{i-1}.data(ori_state_prev_ids)') zeros(3);...
          zeros(3) quat2rotm(x_est{i-1}.data(ori_state_prev_ids)')];
    P_est{i}.data([error_ori_ids error_pos_ids],[error_ori_ids error_pos_ids]) = ...
        Phi * P_est{i-1}.data([error_ori_prev_ids error_pos_prev_ids], [error_ori_prev_ids error_pos_prev_ids]) * Phi' + ...
        G * measurements_config.cov_relative_pose * G';
end

% Plot propagated trajectory
figure('Name','Propagated trajectory'); hold on;
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
grid on
axis equal
for i=1:time.count_time_steps
    ori_state_ids = x_est{i}.sz_idx.ego_idx{i}(x_est{i}.sz_idx.ego.ori_idx);
    pos_state_ids = x_est{i}.sz_idx.ego_idx{i}(x_est{i}.sz_idx.ego.pos_idx);
    p_est = x_est{i}.data(pos_state_ids);
    R_est = quat2rotm(x_est{i}.data(ori_state_ids)');
    scatter3(p_est(1), p_est(2), p_est(3), 15, 'g' , 'filled');
    quiver3(p_est(1), p_est(2), p_est(3),...
            R_est(1,1), R_est(2,1), R_est(3,1),...
            'color', 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver3(p_est(1), p_est(2), p_est(3),...
            R_est(1,2), R_est(2,2), R_est(3,2),...
            'color', 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver3(p_est(1), p_est(2), p_est(3),...
            R_est(1,3), R_est(2,3), R_est(3,3),...
            'color', 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    % Extract and draw ellipse uncertainty
    error_pos_ids = P_est{i}.sz_idx.ego_idx{i}(P_est{i}.sz_idx.ego.pos_idx);
    [U, D] = eigs(P_est{i}.data(error_pos_ids, error_pos_ids))
    k = U(:, 3);
    p = [];
    for tt = 0:0.01:2*pi
        R = cos(tt)*eye(3) + skewsymm(k) * sin(tt) + (1-cos(tt))*k*k';
        v = P_est{i}.data(error_pos_ids, error_pos_ids) * R * U(:,1) + p_est;
        p = [p; v'];
    end
    plot3(p(:,1), p(:,2), p(:,3), 'b');
    grid on
    axis equal
    pause;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Measurement model and update