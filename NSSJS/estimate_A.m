function  [A]    =   estimate_A( D, D0, X, Y, par, Z_ori, sf, sz,A,lables  )
FB=par.fft_B;  
FBC=par.fft_BT;
nl=sz(1);
V1=ConvC(A, FB, nl);
V2=A;
V3=A;
V4=A;
V5=A;

shift=floor(sf/2)-1;
mask = zeros(sz(1), sz(2));
mask(shift+1:sf:sz(1), shift+1:sf:sz(1)) = 1;
maskim = repmat(mask, [1, 1, size(D,2)]);
mask = im2mat(maskim);
D02=D0'*D0;
D2=D'*D;

mu=par.mu;
II_FB = 1./((abs(FB.^2)) +4);
IBD_B= FBC./((abs(FB.^2)) + 4);

G1=zeros(size(V1));
G2=zeros(size(V2));
G3=zeros(size(V3));
G4=zeros(size(V4));
G5=zeros(size(V5));

X = hyperConvert3D(X, sz(1)/sf, sz(2)/sf);
Yhim_up = upsamp_HS(X, sf, sz(1), sz(2), shift);
X2 = im2mat(Yhim_up);
DTX=D'*X2;
D0Y=D0'*Y;
tol=1e-3;
relerr=[];

for i=1:200
%% A
A_k=A;
A = ConvC(V1+G1, IBD_B, nl) + ConvC(V2+G2, II_FB, nl) +  ConvC(V3+G3, II_FB, nl)+ConvC(V4+G4, II_FB, nl)+ConvC(V5+G5, II_FB, nl);
 %% V1
NU1 =ConvC(A, FB, nl)-G1;     
V1_com =(D2+mu*eye(size(D,2)))\(DTX+mu*NU1);               
V1=NU1.*(1-mask)+  V1_com .*mask; 
%% V2

V2=(D02+mu*eye(size(D0,2)))\(D0Y+mu*(A-G2));


%% V3 and V4

for cc=0:1:par.total_patches-1
      index=find(lables==cc);
      ggg=A-G3;
    V3(:,index)= vector_soft_row(ggg(:,index), ...
       par.gamma/mu* (1./( sqrt(sum(ggg(:,index).^2,2)) + eps ))+ eps );
  
   gggg=A - G4;
    V4_temp=gggg(:,index);
    [u,s,v] = svd(V4_temp,'econ');
    ds = diag(s);
    V4(:,index) = u*diag(max( abs(ds) - (par.tau/mu),...
    zeros(size(ds))))*v';
end
V5 = max(A - G5, 0);
%% G1 G2 G3
G1=G1+(V1-ConvC(A, FB, nl));
G2=G2+(V2-A);
G3=G3+(V3-A);
G4=G4+(V4-A);
G5=G5+(V5-A);


relerr(i)=abs(norm(A_k(:)-A(:))/(norm(A_k)+eps));
relerr(1)=1;
    if relerr(i)<tol
        break
    end
a1 =double(im2uint8(D*A)); 
b1 =double(im2uint8(Z_ori)); 
 MSE             =    mean( mean( (a1-b1).^2 ) );
rmse=  sqrt(MSE);

disp( sprintf('Iter %d, RMSE = %3.3f', i, rmse) );  

end





                              

        

     
