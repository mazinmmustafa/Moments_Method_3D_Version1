close all; clear; clc;

data_0d = importdata("mesh/mesh_0d.txt");
data_1d = importdata("mesh/mesh_1d.txt");
data_2d = importdata("mesh/mesh_2d.txt");
data_3d = importdata("mesh/mesh_3d.txt");

[N_0d, ~] = size(data_0d);
[N_1d, ~] = size(data_1d);
[N_2d, ~] = size(data_2d);
[N_3d, ~] = size(data_3d);
% 
figure()
hold on
for i=1:4:N_3d
    plot3([data_3d(i, 1) data_3d(i+1, 1) data_3d(i+2, 1) data_3d(i+3, 1)],...
          [data_3d(i, 2) data_3d(i+1, 2) data_3d(i+2, 2) data_3d(i+3, 2)],...
          [data_3d(i, 3) data_3d(i+1, 3) data_3d(i+2, 3) data_3d(i+3, 3)],...
          "-r", "LineWidth", 1)
end
for i=1:3:N_2d
    plot3([data_2d(i, 1) data_2d(i+1, 1) data_2d(i+2, 1)],...
          [data_2d(i, 2) data_2d(i+1, 2) data_2d(i+2, 2)],...
          [data_2d(i, 3) data_2d(i+1, 3) data_2d(i+2, 3)],...
          "-k", "LineWidth", 1)
end
for i=1:2:N_1d
    plot3([data_1d(i, 1) data_1d(i+1, 1)],...
          [data_1d(i, 2) data_1d(i+1, 2)],...
          [data_1d(i, 3) data_1d(i+1, 3)],...
          "-b", "LineWidth", 1)
end
for i=1:1:N_0d
    plot3([data_0d(i, 1)],...
          [data_0d(i, 2)],...
          [data_0d(i, 3)],...
          "og", "MarkerSize", 5, "MarkerFaceColor", "g")
end
hold off
FontSize = 10;
xlabel("x", "FontSize", FontSize)
ylabel("y", "FontSize", FontSize)
zlabel("z", "FontSize", FontSize)
axis equal
view([30 45])
