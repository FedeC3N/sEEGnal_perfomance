function [pos_elec,size_elec] = draw_head(config)

hold on;
axis equal off

skin_color = [255 243 231]/255;

% Head (big circle)
theta = linspace(0, 2*pi, 100);
r_head = 1; % radius of the head
x_head = r_head * cos(theta);
y_head = r_head * sin(theta);

% Ears (two small circles)
r_ear = 0.2; % radius of the ears
x_ear = r_ear * cos(theta);
y_ear = r_ear * sin(theta);

% Left ear
fill(x_ear - r_head/1.2, y_ear, skin_color); % Slightly darker gray

% Right ear
fill(x_ear + r_head/1.2, y_ear, skin_color);

% Nose (triangle)
x_nose = 2*[-0.07, 0.07, 0]; % X-coordinates of the triangle
y_nose = 2*[0.49, 0.49, 0.55]; % Y-coordinates of the triangle
fill(x_nose, y_nose, skin_color); % Darker gray for nose

% Plot the head at the end
fill(x_head, y_head, skin_color); % Light gray color

% Plot the sensors
% Save memory for the positions and colors
pos_elec = nan(numel(config.complete_channel_labels),2);
size_elec = 50*ones(numel(config.complete_channel_labels),1);

% Read the channels position file
lines = readlines('../shared/head_layouts/elec1005.lay');

% Go through each line
for iline = 1 : size(lines,1)
    
    % Get the channel position in the original "complete_channel_labels"
    % for consistency
    current_line = strsplit(lines(iline),' ');
    
    if size(current_line,2) > 1
        current_channel = current_line(6);
        current_channel_index = ismember(config.complete_channel_labels,current_channel);
        
        % If present, save the position and color
        if sum(current_channel_index) == 1
            
            % Position
            pos_elec(current_channel_index,1) = current_line(2);
            pos_elec(current_channel_index,2) = current_line(3);
            
        end
        
    end
    
end

pos_elec = pos_elec * 2;

scatter(pos_elec(:,1),pos_elec(:,2),size_elec,'MarkerEdgeColor',[0 0 0])
% pos_labels = [pos_elec(:,1)-0.015, pos_elec(:,2)+0.03];
% text(pos_labels(:,1),pos_labels(:,2),config.complete_channel_labels)

end
