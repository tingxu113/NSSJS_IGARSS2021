function [label,rmseee] = fkmeans_vector(X, kkk,iterations,ratio_ergas)
[m n] = size(X);

%%initial centriod.
%pre_centroid = spreadseeds(X, kkk);
rng(10,'twister')
pre_centroid = X(randsample(size(X,1),kkk),:);

%%the SED distance
dist_sam1= bsxfun(@rdivide,pre_centroid *X',sqrt(sum(pre_centroid .^2,2)+eps)*sqrt(sum(X'.^2,1))+eps);
dist_sam2=real(acos(dist_sam1));
dist_sam=dist_sam2*180/pi;
  for i=1:kkk
     E=bsxfun(@minus,X,pre_centroid(i,:));
     E1=sqrt(E.^2);
     E2(:,i)=100*ratio_ergas*sqrt(sum((E1./pre_centroid(i,:)).^2,2))/n;
  end
  ergas=E2';
 % ergas=ergas/max(ergas(:));
  % dist_sam=dist_sam/max(dist_sam(:));
%    for i=1:kkk
%            for j=1:m
%                c=(X(j,:)-pre_centroid (i,:)).^2;
%                cc=sqrt(sum((sqrt(c)./pre_centroid (i,:)).^2));
%                ergas(i,j)=100*ratio_ergas*cc/n;
%            end
%        end
%  dist_rmse=sqrdistance(k,X)';
 % ergas = 100*ratio_ergas*sqrt(((dist_rmse.^2)./ (1/nn*sum(k,2)).^2));   
% distances=dist_sam./sqrt(dist_sam.^2+ergas.^2+eps)+ergas./sqrt(dist_sam.^2+ergas.^2+eps);
  distances=(dist_sam+ergas)/2;
 [~,label] = min(distances);
   centroid=pre_centroid ;
   kk = size(pre_centroid ,1);

%last = 0;

%if ~weight
    % code defactoring for speed
   % while any(label ~= last)
   rmseee=[];
   for iter=1:iterations
       pre_centroid=centroid;
        [~,~,label] = unique(label);
        for mn=1:max(label)
            gg=find(label==mn);
            X_gg=X(gg,:);
          centroid(mn,:)=mean(X_gg);
          % centroid(mn,:)=sqrt(sum(X_gg.^2,1)/size(X_gg,1));
        end     
      %ind = sparse(label,1:n,1,kk,n,n);    
      dist_sam1= bsxfun(@rdivide,centroid*X',sqrt(sum(centroid.^2,2)+eps)*sqrt(sum(X'.^2,1))+eps);
      dist_sam2=real(acos(dist_sam1));
      dist_sam=dist_sam2*180/pi;
      for i=1:kkk
     E=bsxfun(@minus,X,centroid(i,:));
     E1=sqrt(E.^2);
     E2(:,i)=100*ratio_ergas*sqrt(sum((E1./centroid(i,:)).^2,2))/n;
      end
   ergas=E2';        
  % ergas=ergas/max(ergas(:));
 %  dist_sam=dist_sam/max(dist_sam(:));
      distances=(dist_sam+ergas)/2;
      % distances= dist_rmse;
      [~,label] = min(distances);
      label = label';
   rmseee(iter)=norm(pre_centroid-centroid,2)/norm(pre_centroid,2);
    if sum(norm(pre_centroid-centroid,2))<1e-8 %不断迭代直到位置不再变化
    break;
   end
   end

    
  %  dis = ind*(sum(X.^2,2) - 2*max(distances)');

% else
%     while any(label ~= last)
%         % remove empty clusters
%         [~,~,label] = unique(label);
%         % transform label into indicator matrix
%         ind = sparse(label,1:n,weight,k,n,n);
%         % compute centroid of each cluster
%         centroid = (spdiags(1./sum(ind,2),0,k,k)*ind)*X;
%         % compute distance of every point to each centroid
%         dist_sam1= bsxfun(@rdivide,centroid*X',sqrt(sum(centroid.^2,2)+eps)*sqrt(sum(X'.^2,1))+eps);
%        %tmp_sam=bsxfun(@rdivide,k*X',sqrt(sum(k.^2,2)+eps)*sqrt(sum(X'.^2,1))+eps);
%         dist_sam2=real(acos(dist_sam1));
%         dist_sam=dist_sam2*180/pi;
%      % dist_rmse=norm(data(i,:)-centroids(j,:),2);
%         dist_rmse=sqrdistance(centroid,X)';
%         distances=dist_sam./sqrt(dist_sam.^2+dist_rmse.^2+eps)+dist_rmse./sqrt(dist_sam.^2+dist_rmse.^2+eps);
%        % distances = bsxfun(@minus,centroid*X',0.5*sum(centroid.^2,2));
%         % assign points to their nearest centroid
%         last = label;
%         [~,label] = max(distances);
%         label = label';
%     end
%     dis = ind*(sum(X.^2,2) - 2*max(distances)');
% end
% label = label';

% Code below this line reused from the file exchange submission K-means++
% (http://www.mathworks.com/matlabcentral/fileexchange/28901-k-means) in
% accordance with the license:
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
% Copyright (c) 2010, Michael Chen
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% 
function D = sqrdistance(A, B)
% Square Euclidean distances between all sample pairs
% A:  n1 x d data matrix
% B:  n2 x d data matrix
% WB: n2 x 1 weights for matrix B
% D: n2 x n1 pairwise square distance matrix
%    D(i,j) is the squared distance between A(i,:) and B(j,:)
% Written by Michael Chen (sth4nth@gmail.com). July 2009.
n1 = size(A,1); n2 = size(B,2);
m = (sum(A,1)+sum(B,1))/(n1+n2);
A = bsxfun(@minus,A,m);
B = bsxfun(@minus,B,m);
D = full((-2)*(A*B'));
D = bsxfun(@plus,D,full(sum(B.^2,2))');
D = bsxfun(@plus,D,full(sum(A.^2,2)))';
end

function [S, idx] = spreadseeds(X, k)
% X: n x d data matrix
% k: number of seeds
% reference: k-means++: the advantages of careful seeding.
% by David Arthur and Sergei Vassilvitskii
% Adapted from softseeds written by Mo Chen (mochen@ie.cuhk.edu.hk), 
% March 2009.
[n,d] = size(X);
idx = zeros(k,1);
S = zeros(k,d);
D = inf(n,1);
idx(1) = ceil(n.*rand);
S(1,:) = X(idx(1),:);
for i = 2:k
    D = min(D,sqrdistance(S(i-1,:),X));
    idx(i) = find(cumsum(D)/sum(D)>rand,1);
    S(i,:) = X(idx(i),:);
end
end

end
