clc;clear all;close all;

M = 16; 
G = 10;
c = 3e8; 
f = 5.8e9;   
lambda = c/f; 
d = lambda/2; 
angles = [-60,-20,10,15,40]; 
t = 0:1e-11:10000*1e-11;
S = [exp(-1j*2*pi*1e9*t);exp(-1j*2*pi*2e9*t);exp(-1j*2*pi*3e9*t);...
    exp(-1j*2*pi*4e9*t);exp(-1j*2*pi*5e9*t)];  %非相干
S0 = [exp(-1j*2*pi*1e9*t);exp(-1j*2*pi*2e9*t);exp(-1j*2*pi*2e9*t);...
    exp(-1j*2*pi*2e9*t);exp(-1j*2*pi*3e9*t)];  %相干
A = exp(-1j*2*pi*(0:M-1)'*d/lambda*sind(angles)); 
X = A*S0;
X = awgn(X, 10, 'measured');
Rx = X*X';

theta=-90:0.05:90; 
a=exp(-1j*2*pi*(0:M-1)'*d/lambda*sind(theta));  % MUSIC
aF=exp(-1j*2*pi*(0:G-1)'*d/lambda*sind(theta)); % FSS MUSIC

figure;
angles_str = strjoin(arrayfun(@(x) [num2str(x), '°'], angles, 'UniformOutput', false), ', ');
sgtitle(['信号2、3、4相干，Angles = ', angles_str]);
% CBF
subplot(4,1,1);
w1 = a;
p1 = abs(w1'*Rx*w1);
plot(theta,10*log10(p1/max(p1)),'linewidth',2);
title('CBF');xlabel('入射角度(°)');ylabel('功率(dB)');
hold on;
for i = 1:length(angles)
    x_line = angles(i);
    line([x_line x_line], ylim, 'Color', 'r', 'LineStyle', '--');  % 绘制红色的垂直线
end
hold off;

% MVDR
subplot(4,1,2)
p2 = abs(diag(1./(a'/Rx*a)));
plot(theta,10*log10(p2/max(p2)),'linewidth',2);
title('MVDR');xlabel('入射角度(°)');ylabel('功率(dB)');
hold on;
for i = 1:length(angles)
    x_line = angles(i);
    line([x_line x_line], ylim, 'Color', 'r', 'LineStyle', '--');  % 绘制红色的垂直线
end
hold off;

% MUSIC
subplot(4,1,3)
[U, L] = eig(Rx);
UU = fliplr(U);
Un = UU(:,(length(angles)+1):M);
p3 = abs(diag(1./(a'*(Un*Un')*a)));
plot(theta,10*log10(p3/max(p3)),'linewidth',2);
title('MUSIC');xlabel('入射角度(°)');ylabel('功率(dB)');
hold on;
for i = 1:length(angles)
    x_line = angles(i);
    line([x_line x_line], ylim, 'Color', 'r', 'LineStyle', '--');  % 绘制红色的垂直线
end
hold off;

% FSS MUSIC
RxF = zeros(G);
for i = 1:M-G-1
    RxF = RxF+Rx(i:i+G-1,i:i+G-1);
end
RxF = RxF/(M-G-1);
subplot(4,1,4)
[UF, LF] = eig(RxF);
UUF = fliplr(UF);
UnF = UUF(:,(length(angles)+1):G);
p4 = abs(diag(1./(aF'*(UnF*UnF')*aF)));
plot(theta,10*log10(p4/max(p4)),'linewidth',2);
title('FSS MUSIC');xlabel('入射角度(°)');ylabel('功率(dB)');
hold on;
for i = 1:length(angles)
    x_line = angles(i);
    line([x_line x_line], ylim, 'Color', 'r', 'LineStyle', '--');  % 绘制红色的垂直线
end
hold off;
