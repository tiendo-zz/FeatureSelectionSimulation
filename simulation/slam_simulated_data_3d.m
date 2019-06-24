%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load configuration settings
% TODO: Move implementation of trajectory generators in seperate file

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
measurements_config.cov_relative_pose = [0.02^2, 0, 0,     0, 0, 0;
                                         0, 0.02^2, 0,     0, 0, 0;
                                         0, 0, 0.02^2,     0, 0, 0;
                                         0, 0, 0,      0.1^2, 0, 0;
                                         0, 0, 0,      0, 0.1^2, 0;
                                         0, 0, 0,      0, 0, 0.1^2];                                         
                                                                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate and visualize data for a slam experiment

% Real world simulation: Generate data for a slam experiment
[time, trajectory, measurements] = real_world_simulation(time_config, trajectory_config, measurements_config);

% Visualize the trajectory 
figure('Name','EKF animation'); hold on;
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





% Initialize the animation figure
if animation
    animation_fig = figure('Name','EKF animation');
    disp(['Time step: ', num2str(0.0), ' Cloned pose: ',num2str(1)]);
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    grid on
    axis equal
    view(-35,45);
    scatter3(G_p_S_true(1,animation_indices(1)),...
             G_p_S_true(2,animation_indices(1)),...
             G_p_S_true(3,animation_indices(1)), 15, 'g' , 'filled');    
    for anim_idx = 2:length(animation_indices)
        idx_end = animation_indices(anim_idx);
        idx_start = animation_indices(anim_idx-1);
        figure(animation_fig); hold on;
        scatter3(G_p_S_true(1,idx_end), G_p_S_true(2,idx_end), G_p_S_true(3,idx_end), 15, 'g' , 'filled');
        plot3(G_p_S_true(1,idx_start:idx_end),...
              G_p_S_true(2,idx_start:idx_end),...
              G_p_S_true(3,idx_start:idx_end), 'g'); 
        axis equal
    end
    scatter3(x_est(pos_idx(1)), x_est(pos_idx(2)), x_est(pos_idx(3)), 15, 'b' , 'filled');
    R_est = quat2rotm(x_est(ori_idx)');
    quiver3(x_est(pos_idx(1)), x_est(pos_idx(2)), x_est(pos_idx(3)),...
            R_est(1,1), R_est(2,1), R_est(3,1),...
            'color', 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver3(x_est(pos_idx(1)), x_est(pos_idx(2)), x_est(pos_idx(3)),...
            R_est(1,2), R_est(2,2), R_est(3,2),...
            'color', 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver3(x_est(pos_idx(1)), x_est(pos_idx(2)), x_est(pos_idx(3)),...
            R_est(1,3), R_est(2,3), R_est(3,3),...
            'color', 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
end


figure(animation_fig); hold on;
scatter3(x_est_rt(pos_idx(1),idx_clone),...
                 x_est_rt(pos_idx(2),idx_clone),...
                 x_est_rt(pos_idx(3),idx_clone), 15, 'b' , 'filled');
        R_est = quat2rotm(x_est_rt(ori_idx,idx_clone)');
        quiver3(x_est_rt(pos_idx(1),idx_clone),...
                x_est_rt(pos_idx(2),idx_clone),...
                x_est_rt(pos_idx(3),idx_clone),...
                R_est(1,1), R_est(2,1), R_est(3,1),...
                'color', 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
        quiver3(x_est_rt(pos_idx(1),idx_clone),...
                x_est_rt(pos_idx(2),idx_clone),...
                x_est_rt(pos_idx(3),idx_clone),...
                R_est(1,2), R_est(2,2), R_est(3,2),...
                'color', 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
        quiver3(x_est_rt(pos_idx(1),idx_clone),...
                x_est_rt(pos_idx(2),idx_clone),...
                x_est_rt(pos_idx(3),idx_clone),...
                R_est(1,3), R_est(2,3), R_est(3,3),...
                'color', 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
        