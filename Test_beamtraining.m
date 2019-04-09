clear;clc;
w_r = (-pi/2 : 0.01 : pi/2);
sin_w_r = [-1:0.0001:1];
w_r = asin(sin_w_r);
N = 32;  %天线数量
N_RF = 4;
N_block = N / N_RF;
N_beam = N_block * (N_RF-2);
CodeBook = dftmtx(N);
% subCodeBook = dftmtx(N/4);
select_en = 1;  % 1:not rand ; 0:rand

%  %% DFT codebook
% A_DFT = dftmtx(N);
% 
% %% RF部分
% F_RF_nperm = A_DFT;
% index_F = randperm(N);    % random permutation
% F_RF_perm = F_RF_nperm(:,index_F);
% 
% %% BB部分
% F_BB_q = dftmtx(N_RF);
% F_BB_q = F_BB_q(:,1:N_beam/N_block);
% F_BB = kron(eye(N_block),F_BB_q);
% F = F_RF_perm*F_BB;
% F = sqrt(N_beam)/norm(F,'fro')*F;

% 生成扫描矢量
%base_vector = exp(1i*pi*sin(w_r)'*(0:N-1));
base_vector = exp(1i*pi*sin_w_r'*(0:N-1));
if select_en==1
%     angle = -pi/6;  %波束对齐角度
%     target_vector = exp(1i*pi*sin(angle)*[0:1:N-1]); 

%     target_vector = [subCodeBook(:,1).' subCodeBook(:,2).' subCodeBook(:,3).' subCodeBook(:,4).'];

   % target_vector = CodeBook(:,2).';(CodeBook(:,1).'+ CodeBook(:,2).'+ CodeBook(:,3).'+ CodeBook(:,4).'); sum(CodeBook.');
   % target_vector =(CodeBook(:,1).'+ CodeBook(:,2).'+ CodeBook(:,3).'); sum(CodeBook.');
     target_vector =sum(CodeBook(:,[4 4]).');


elseif select_en==2
    target_vector = exp(1i*2*pi*rand(1,N));
end
%target_vector([8:24]) = 0;
f_r = abs(target_vector*base_vector')/N;
figure(3);plot(abs(f_r)');
figure(1);
polarplot(w_r+pi/2, f_r, '-')   %加pi/2，将图像移动到上半圆
ax = gca;
ax.ThetaLim = [0 180];
ax.ThetaTickLabel = {'-90','-60','-30','0','30','60','90'};
hold on;

  
