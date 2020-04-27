clear all;
close all;
clc;

data = csvread('data.csv', 1, 0);
datatable = readtable('data.csv', 'ReadVariableNames', false);
numRows = height(datatable);
f1 = figure;
for N = 1 : numRows-1
    if mod(N, 2) == 0
        clf;
        plot(0, 0, "r*");
        hold on;
        x = data(N, 2:11);
        y = data(N+1,2:11);
        plot(x, y, "bo");
        hold on;
        xs = data(N, 1);
        ys = data(N+1, 1);
        plot(xs, ys, "gp");
        axis([-600 600 -600 600])
        title("System1: Satellite launch T+" + N/2 + " days", 'FontSize', 10);
        xlabel('X', 'FontSize', 10);
        ylabel('Y', 'FontSize', 10);
        hold on;
        grid on;
        pause(.005);
    end
end

data2 = csvread('data2.csv', 1, 0);
datatable2 = readtable('data2.csv', 'ReadVariableNames', false);
numRows2 = height(datatable2);
f2 = figure;
for N2 = 1 : numRows2-1
    if mod(N2, 2) == 0
        clf;
        plot(0, 0, "r*");
        hold on;
        x2 = data(N2, 2:5);
        y2 = data(N2+1,2:5);
        plot(x2, y2, "bo");
        hold on;
        xs2 = data2(N2, 1);
        ys2 = data2(N2+1, 1);
        plot(xs2, ys2, "gp");
        axis([-300 300 -300 300])
        title("System2: Satellite launch T+" + N2/2 + " days", 'FontSize', 10);
        xlabel('X', 'FontSize', 10);
        ylabel('Y', 'FontSize', 10);
        hold on;
        grid on;
        pause(.005);
    end
end
