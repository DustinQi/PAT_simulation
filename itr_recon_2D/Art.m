function [X]=Art(W,b,X,Niter,e)
% =============================================================
% ==== J. Ripoll's codes for light diffusion ==================
% =============================================================
% [X,err,br]=art(W,b,X,Niter,e);
% -----------------------------------------
% Inverts matrix X by using the Algebraic Reconstruction Technique
% (ART) - See Kak & Slaney, Ch. 7
% W*X = b
% | W11 W12 W13 ... W1N| |X1|   |b1|
% | W21 W22 W23 ... W2N| |X2|   |b2|
% | W31 W32 W33 ... W3N| |X3| = |b3|
% | :    :   :  ...  : | |: |   |: |
% | WM1 WM2 WM3 ... WMN| |XM|   |bN|
% -------------------------------------------------------------
% Input:
% -----
% W......| The weight matrix.
% b......| The measurement matrix.
% X......| The matrix to invert (with born/rytov its the perturbation)
% Niter..| Number of iterations
% e......| Relaxation factor.
%
% Output:
% -------
% X......| Reconstructed matrix.
% err....| Error commited in the reconstruction
% br.....| Final solution to W*Xrec.
% -------------------------------------------------------------
%x,b输入为行向量W*X'=b'
[N,M]=size(W);
%Norm = diag(W*W');
Norm = zeros(N,1);
for ii=1:N
   Norm(ii) = sum( W(ii,:).^2 ); 
end

for ii=1:Niter,
	for jj=1:N,
      	Q = (X*W(jj,:)');
		Delta = (b(jj)-Q)./Norm(jj)*W(jj,:);
		X = X + e*Delta;          
      end;
      %X = X.*(1+sign(X))/2.*abs(sign(X));
%       X(X<0) = 0;
      if(mod(ii,10)==0)
           % echo(ii)
           ii
      end
end;
     
