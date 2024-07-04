close all; clear; clc;
%%

RCS1 = load("../data/RCS1.txt");

figure()
hold on
plot(RCS1(:, 1), 10*log10(RCS1(:, 2)))
hold off

RCS2 = load("../data/RCS2.txt");

figure()
hold on
plot(RCS2(:, 1), 10*log10(RCS2(:, 3)))
hold off