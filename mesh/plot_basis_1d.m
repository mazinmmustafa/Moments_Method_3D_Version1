close all; clear; clc;
%
data_1d = load("basis/basis_1d.txt");
%
[N_1d, ~] = size(data_1d);
% 
figure()
FontSize = 10;
xlabel("x", "FontSize", FontSize)
ylabel("y", "FontSize", FontSize)
zlabel("z", "FontSize", FontSize)
axis equal
view([30 45])
hold on
for i=1:N_1d
    plot3([data_1d(i, 1) data_1d(i, 1+3)],...
          [data_1d(i, 2) data_1d(i, 2+3)],...
          [data_1d(i, 3) data_1d(i, 3+3)],...
          "-r", "LineWidth", 1)
    plot3([data_1d(i, 1+3+3) data_1d(i, 1+3)],...
          [data_1d(i, 2+3+3) data_1d(i, 2+3)],...
          [data_1d(i, 3+3+3) data_1d(i, 3+3)],...
          "-b", "LineWidth", 1)
      input(" ");
end
hold off

