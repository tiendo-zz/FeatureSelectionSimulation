function [time, trajectory, world, measurements] = real_world_simulation(time_config, trajectory_config, world_config, camera_config, measurements_config)
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


% True 3D features' positions w.r.t the first camera as inverse depth
world.point_feat_Wposition = cell(world_config.point_feature_num,1);
world.point_feat_invd = cell(world_config.point_feature_num,1);
k = 1;
while( k <= world_config.point_feature_num )
    tt = 2*pi*rand();
    phi = pi*rand() - pi/2;
    W_p_f = (0.5 + 0.5*rand())*world_config.radius*[cos(phi)*cos(tt);cos(phi)*sin(tt);sin(phi)]; % <--- w.r.t the world    
    flag = 0;    
    observation_count = 0;
    for i = 1:time.count_time_steps
        W_p_Ci = trajectory_config.position(time.absolute_time_true(i));
        W_R_Ci = trajectory_config.orientation(time.absolute_time_true(i));
        Ci_p_f = W_R_Ci'*(W_p_f - W_p_Ci);
        fov = atan2(norm(Ci_p_f(1:2)), Ci_p_f(3));
        if fov < camera_config.fov/2
            if flag == 0
                flag = 1;
                First_view_cam_idx = i;                
                world.point_feat_position{k} = [Ci_p_f;First_view_cam_idx];
                world.point_feat_Wposition{k} = W_p_f;
                
                world.point_feat_invd{k} = zeros(4,1);
                world.point_feat_invd{k}(3) = 1/norm(Ci_p_f(1:3));
                world.point_feat_invd{k}(1) = atan2(Ci_p_f(2),Ci_p_f(1));
                world.point_feat_invd{k}(2) = atan2(Ci_p_f(3),(Ci_p_f(1)^2+Ci_p_f(2)^2)^0.5);
                world.point_feat_invd{k}(4) = First_view_cam_idx;
            end            
            observation_count = observation_count + 1;
        end
    end
    if observation_count == time.count_time_steps
        k = k + 1; % Next feature generating
    end
end


% Pixel projection measurement
% Create feature projections as measurement
measurements.corner_pixel_true = cell(time.count_time_steps, 1);
measurements.corner_pixel = cell(time.count_time_steps, 1);
measurements.corner_homogeneous = cell(time.count_time_steps, 1);

for i = 1:time.count_time_steps    
    W_p_Ci = trajectory_config.position(time.absolute_time_true(i));
    W_R_Ci = trajectory_config.orientation(time.absolute_time_true(i));
    % TODO: add visibility graph and only project features that are visible
    % on camera i
    measurements.corner_pixel{i} = zeros(3, world_config.point_feature_num);
    measurements.corner_homogeneous{i} = zeros(4, world_config.point_feature_num);
    for l = 1:world_config.point_feature_num
        Ci_p_f = W_R_Ci'*(world.point_feat_Wposition{l} - W_p_Ci);
        
        % [u v 1 feat_id], feat_id needs to be fixed corresponds to
        % visibility graph
        measurements.corner_homogeneous{i}(1:3, l) = Ci_p_f./Ci_p_f(3);
        measurements.corner_homogeneous{i}(4, l) = l;
        % [x y feat_id], feat_id needs to be fixed corresponds to
        % visibility graph
        measurements.corner_pixel_true{i}(1:2, l) = camera_config.fc*measurements.corner_homogeneous{i}(1:2, l)+...
                                                    camera_config.cc;
        measurements.corner_pixel_true{i}(3, l) = l;
        measurements.corner_pixel{i}(:, l) = measurements.corner_pixel_true{i}(:, l);
        measurements.corner_pixel{i}(1:2, l) = measurements.corner_pixel_true{i}(1:2, l) + randn(2,1)*camera_config.sigma_pixel;
    end    
end
                                                 
% Perturb features' positions


% Visual measurement graph


end
