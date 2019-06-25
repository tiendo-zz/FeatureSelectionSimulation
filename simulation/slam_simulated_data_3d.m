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
measurements_config.cov_relative_pose = [0.2^2, 0, 0,     0, 0, 0;
                                         0, 0.2^2, 0,     0, 0, 0;
                                         0, 0, 0.2^2,     0, 0, 0;
                                         0, 0, 0,      1^2, 0, 0;
                                         0, 0, 0,      0, 1^2, 0;
                                         0, 0, 0,      0, 0, 1^2];   
                                     
delay_line = 10;
                                                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate and plot/visualize data for a slam experiment

% Real world simulation: Generate data for a slam experiment
[time, trajectory, measurements] = real_world_simulation(time_config, trajectory_config, measurements_config);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate and plot/visualize data for a slam experiment

state.ori_size = 4; 
state.ori_idx = 1:4;
state.pos_size = 3; 
state.pos_idx = 5:7;
state.pose_size = state.ori_size + state.pos_size;
state.pose_idx = [state.ori_idx, state.pos_idx];
state.pose_count = 1;

error_state.ori_size = 3; 
error_state.ori_idx = 1:3;
error_state.pos_size = 3; 
error_state.pos_idx = 4:6;
error_state.pose_size = error_state.ori_size + error_state.pos_size;
error_state.pose_idx = [error_state.ori_idx, error_state.pos_idx];
error_state.pose_count = 1;

x_est = cell(1, time.count_time_steps);
P_est = cell(1, time.count_time_steps);
% estimated state and covariance
for i= 1:time.count_time_steps
    window_size = i - max(i-delay_line,1) + 1;
    
    x_est{i}.data = zeros(window_size*state.pose_size,1);
    x_est{i}.sz_idx = state;
    x_est{i}.sz_idx.count_states = window_size;
    % Fill the pose index
    for j=1:window_size
        
    end
    
    P_est{i}.data = zeros(window_size*error_state.pose_size,window_size*error_state.pose_size);
    P_est{i}.sz_idx = error_state;
    P_est{i}.count_states = window_size;
    % Fill the pose index
    for j=1:window_size
        
    end
end




