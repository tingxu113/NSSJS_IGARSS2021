function    [HSI_res1]  =   NSSJS_HSI_SR( Z_ori, sf, par,sz,X,Y,D)
 Y1=normalize(Y');
 B=sparsepca(Y1);
 Y1=Y1*B;
 Y1=reshape(Y1,sz(1),sz(2));

[labels,~]=suppixel(Y1,par.total_patches);
A         =   zeros( par.K, size(Y,2) );

for i=1
 [A]   =   estimate_A(D, par.P*D, X, Y, par, Z_ori, sf, sz,A,labels );
%  a1 =double(im2uint8(D*A)); 
% b1 =double(im2uint8(Z_ori)); 
%  MSE2             =    mean( mean( (a1-b1).^2 ) );
% rmse1(i)=  sqrt(MSE2)
% a1 =double(im2uint8(D*A)); 
% b1 =double(im2uint8(Z_ori)); 
% MSE2 = mean( mean( (a1-b1).^2 ) );
% rmse(i)=  sqrt(MSE2)
end
HSI_res1=D*A;



