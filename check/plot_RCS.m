close all; clear; clc; FontSize=24;
%%

RCS1 = importdata("../data/RCS1.txt");

figure()
hold on
plot(RCS1(:, 1)*180/pi, 10*log10(RCS1(:, 2)), "-k", "LineWidth", 1)
hold off
xlim([0 180])
pbaspect([1.4 1 1])

RCS2 = importdata("../data/RCS2.txt");

figure()
hold on
plot(RCS2(:, 1)*180/pi, 10*log10(RCS2(:, 3)), "-k", "LineWidth", 1)
hold off
xlim([0 180])
pbaspect([1.4 1 1])