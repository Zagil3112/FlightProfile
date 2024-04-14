%% Orbit Graph
image_file = 'land_ocean_ice_2048.jpg';
cdata = imread(image_file);
figure(1);
hold on
[x,y,z] = ellipsoid(0, 0, 0, RE, RE, RE, 30); % Centre of Earth
globe = surface(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]); % Plot Globe
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 0.9, 'EdgeColor', 'none'); % Globe Characteristics
grid on;
axis equal;
for i = 1:3:3*xsize
    plot3(Rt(:,i),Rt(:,i+1),Rt(:,i+2),'.');
end