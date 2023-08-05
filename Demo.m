%--------------------------------------------------------------------------
% DEMO: A VARIATIONAL APPROACH WITH NONLOCAL SELF-SIMILARITY AND
% JOINT-SPARSITY FOR HYPERSPECTRAL IMAGE SUPER-RESOLUTION
% 
% References:
% @inproceedings{xu2021variational,
%   title={A Variational Approach with Nonlocal Self-Similarity and Joint-Sparsity for Hyperspectral Image Super-Resolution},
%   author={Xu, Ting and Huang, Ting-Zhu and Chen, Yong and Huang, Jie and Deng, Liang-Jian},
%   booktitle={2021 IEEE International Geoscience and Remote Sensing Symposium IGARSS},
%   pages={2444-2447},
%   year={2021},
%   organization={IEEE}

% Author: Ting Xu
% Email : tingxu113@163.com
%--------------------------------------------------------------------------
clear
clc

addpath('functions', 'NSSJS', 'superpixel')
S=imread('original_rosis.tif');
%% HR-HSI
S=double(S);
S=S(1:256,1:256,11:end);
S=S/max(S(:));
[M,N,L] = size(S);
sz=[M N];
%% the downsampling matrix of the  spectral mode
 F=load('R.mat');
 F=F.R;
 T=F(:,11:end);
 for band = 1:size(T,1)
     div = sum(T(band,:));
     for i = 1:size(T,2)
         T(band,i) = T(band,i)/div;
     end
 end
%%  simulate LR-HSI 
S1 = hyperConvert2D(S);
sf=16;
s0=sf/2;
BW=ones(16,1)/16;
BW1=psf2otf(BW,[M 1]);
S_w=ifft(fft(S).*repmat(BW1,1,N,L)); 

BH=ones(16,1)/16;
BH1=psf2otf(BH,[N 1]);
aa=fft(permute(S_w,[2 1 3]));
S_h=(aa.*repmat(BH1,1,M,L));
S_h= permute(ifft(S_h),[2 1 3]);  
 
Y_h=S_h(s0:sf:end,s0:sf:end,:);
X=hyperConvert2D(Y_h); 
%% simulate HR-MSI
Y = T*S1;
%% NSSJS
 kernel_type  =  {'uniform_blur'};
 par = Parameters_setting(sf, kernel_type, sz);
 par.total_patches =160;
 par.K=150;
 par.P=T;
 par.gamma=1e-6;
 par.tau=1e-4;
 par.mu=1e-3;
 tic
  [D,~]  =  Nonnegative_DL(X, par ); 
  [M1] = NSSJS_HSI_SR(S1,sf,par,sz,X,Y,D);
 RUN_time = toc;
 Z5=hyperConvert3D(M1, M, N);
 [psnr5,rmse5, ergas5, sam5, uiqi5,ssim5,DD5] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z5)), 0, 1.0/sf);








