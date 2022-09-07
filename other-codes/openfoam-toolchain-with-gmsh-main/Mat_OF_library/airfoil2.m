% function [xx, yy, x, y] = airfoil2(ss, yTE, x1, ch1, th1, tw1, x2, ch2, th2, tw2)

ss  = 1;
yTE = -0.01;
x1  = 
ch1 =
th1 = 
tw1 = 
x2  =
ch2 = 
th2 =
tw2 =


positions = zeros(5, 2);
positions(1, 1) = 1;                positions(1, 2) = yTE(1) / 2;
positions(2, 1) = x2(1) + tw2 / 2;  positions(2, 2) = ch2(1) + th2(1) / 2;
positions(3, 1) = x1(1) + tw1 / 2;  positions(3, 2) = ch1(1) + th1(1) / 2;
positions(4, 1) = 0;                positions(4, 2) = 0;
positions(5, 1) = x1(1) - tw1 / 2;  positions(5, 2) = ch1(1) - th1(1) / 2;
positions(6, 1) = x2(1) - tw2 / 2;  positions(6, 2) = ch2(1) - th2(1) / 2;
positions(7, 1) = 1;                positions(7, 2) = -yTE(1) / 2;

x = positions(:, 1);
y = positions(:, 2);
s = 0 : 1 / (length(positions) - 1) : 1;

%ss = 0 : 0.01 : 1;
xx = spline(s, x, ss);
yy = spline(s, y, ss);