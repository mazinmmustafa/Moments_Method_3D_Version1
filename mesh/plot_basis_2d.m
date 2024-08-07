close all; clear; clc;
%
data_2d = load("basis/basis_2d.txt");
%
[N_2d, ~] = size(data_2d);
% 
figure()
FontSize = 10;
xlabel("x", "FontSize", FontSize)
ylabel("y", "FontSize", FontSize)
zlabel("z", "FontSize", FontSize)
axis equal
view([30 45])
hold on
for i=1:N_2d
      if 1
    plot3([data_2d(i, 1) data_2d(i, 1+3) data_2d(i, 1+3+3) data_2d(i, 1+3+3+3)],...
          [data_2d(i, 2) data_2d(i, 2+3) data_2d(i, 2+3+3) data_2d(i, 2+3+3+3)],...
          [data_2d(i, 3) data_2d(i, 3+3) data_2d(i, 3+3+3) data_2d(i, 3+3+3+3)],...
          "-k", "LineWidth", 1)
    plot3([data_2d(i, 1) data_2d(i, 1+3+3)],...
          [data_2d(i, 2) data_2d(i, 2+3+3)],...
          [data_2d(i, 3) data_2d(i, 3+3+3)],...
          "-k", "LineWidth", 1)
    plot3([data_2d(i, 1+3+3+3) data_2d(i, 1+3)],...
          [data_2d(i, 2+3+3+3) data_2d(i, 2+3)],...
          [data_2d(i, 3+3+3+3) data_2d(i, 3+3)],...
          "-k", "LineWidth", 1)
    plot3([data_2d(i, 1) (data_2d(i, 1+3)+data_2d(i, 1+3+3))/2],...
          [data_2d(i, 2) (data_2d(i, 2+3)+data_2d(i, 2+3+3))/2],...
          [data_2d(i, 3) (data_2d(i, 3+3)+data_2d(i, 3+3+3))/2],...
          "-r", "LineWidth", 1)
    plot3([data_2d(i, 1+3+3+3) (data_2d(i, 1+3)+data_2d(i, 1+3+3))/2],...
          [data_2d(i, 2+3+3+3) (data_2d(i, 2+3)+data_2d(i, 2+3+3))/2],...
          [data_2d(i, 3+3+3+3) (data_2d(i, 3+3)+data_2d(i, 3+3+3))/2],...
          "-b", "LineWidth", 1)
%     input(" ");
      end
end
hold off

