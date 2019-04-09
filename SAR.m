clc;
clear;
close all;
c = 3e8;%光速
tp = 0.2e-6;%脉冲宽度
B = 150e6;%发射信号带宽
fs = 200e6;%采样率
kr = B/tp;%调频率
H = 3000;%高度
lambda = 0.05;%载波波长
fc = c/lambda;%载波频率
A = 0.025;%方位向波束宽度
Rmin = 15000;%初始采样距离
PRT = 1/1000;%脉冲重复时间
v0 = 100;%飞机飞行速度

% 目标参数
dis = [15200, 15200, 15250];
xt = [0, 30, 0];
Ls = A*dis;
% dis = [15200];
% xt = [0];


Ns = 5000;%方位向采样次数
t = -tp*10:1/fs:tp*10;
Nf = length(t);%距离向采样点数
x0 = -Ns*PRT*v0/2;%初始采样开始位置
tf = t+Rmin*2/c;
n = 1:Ns;
im = zeros(Ns, Nf);

% 产生回波信号
for i = 1:length(xt)
Rn = sqrt((x0+n*PRT*v0-xt(i)).^2+dis(i)^2);
tau = 2*Rn/c;
t_ = ones(Ns,1)*tf-tau'*ones(1,Nf);
phase = pi*kr*t_.^2-(4*pi/lambda*Rn')*ones(1,Nf);
im = im+exp(1i*phase).*rectpuls(t_/tp).*(rectpuls((x0+n*PRT*v0-xt(i))/Ls(i))'*ones(1,Nf));
end
image(mat2gray(angle(im))*255);
colormap gray;
num = 10;
xl = linspace(Rmin-5*tp*c, Rmin+5*tp*c, num)/10^3;
set(gca,'XTick',linspace(1, Nf, num));
set(gca,'XTickLabel',regexp(num2str(xl, '%.3f\n'), '\n', 'split'));
yl = linspace(x0, x0+Ns*PRT*v0, num);
set(gca,'YTick',linspace(1, Ns, num));
set(gca,'YTickLabel',regexp(num2str(yl, '%.1f\n'), '\n', 'split'));
title('三点目标回波信号');

% 距离向压缩
h_ref = exp(1i*pi*kr*t.^2).*rectpuls(t/tp);
h_ref = h_ref.* (hamming(Nf).');
H_ref = fftshift(fft(fftshift((ones(Ns, 1) * h_ref).'))).';
im1 = fftshift(fft(fftshift(im.'))).'.* conj(H_ref);
im2 = fftshift(ifft(fftshift(im1.'))).';
figure;
image(mat2gray(abs(im2)) * 255);
colormap gray;
set(gca,'XTick',linspace(1, Nf, num));
set(gca,'XTickLabel',regexp(num2str(xl, '%.3f\n'), '\n', 'split'));
set(gca,'YTick',linspace(1, Ns, num));
set(gca,'YTickLabel',regexp(num2str(yl, '%.1f\n'), '\n', 'split'));
title('三点目标距离向压缩');

% 徙动校正
Ka = -2*v0^2/(lambda*Rmin);
fr = (-Nf/2:(Nf/2 - 1))*fs/Nf;
ff = ones(Ns,1)*fr;
fda = (-Ns/2:(Ns/2-1))/PRT/Ns;
fu = fda'*ones(1,Nf);
D_ref = exp(1i*pi/fc^2/Ka*(fu.*ff).^2-1i*pi/fc/Ka*fu.^2.*ff);
im3 = fftshift(fft(fftshift(im1)));
im3 = im3.*D_ref;
im3 = fftshift(ifft(fftshift(im3.'))).';
im3 = fftshift(ifft(fftshift(im3)));
figure;
image(mat2gray(abs(im3))*255);
colormap gray;
set(gca,'XTick',linspace(1,Nf,num));
set(gca,'XTickLabel',regexp(num2str(xl,'%.3f\n'),'\n','split'));
set(gca,'YTick',linspace(1,Ns,num));
set(gca,'YTickLabel',regexp(num2str(yl,'%.1f\n'),'\n','split'));
title('三点目标距离徙动校正');

% 方位向压缩
tr = (-Ns/2:(Ns/2-1))*PRT;
a_ref = exp(1i*pi*Ka*tr.^2).*rectpuls(tr/A/Rmin*v0);
A_ref = fftshift(fft(fftshift(a_ref)));
A_ref = conj(A_ref).'* ones(1,Nf);
im4 = fftshift(fft(fftshift(im3))).*A_ref;
im4 = fftshift(ifft(fftshift(im4)));
figure;
image(mat2gray(abs(im4))*255);
colormap gray;
set(gca,'XTick',linspace(1,Nf,num));
set(gca,'XTickLabel',regexp(num2str(xl,'%.3f\n'),'\n','split'));
set(gca,'YTick',linspace(1,Ns,num));
set(gca,'YTickLabel',regexp(num2str(yl,'%.1f\n'),'\n','split'));
title('三点目标方位向压缩');
