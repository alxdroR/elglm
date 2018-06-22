function t = toeplitzblk(c,r)
%TOEPLITZBLK Block Toeplitz matrix.
%   TOEPLITZBLK(C,R) is a non-symmetric Toeplitz matrix having C as its
%   first column matrix and R as its first row matrix.   
%
%   Class support for inputs C,R:
%      float: double, single
%

% adr 
% revised: 8/10/2013

% number of rows in block
[br,M] = size(r);
% number of columns
[T,bc] = size(c);

if r(:,1:bc) ~= c(1:br,:)
    warning(message('toeplitzblk:DiagonalConflict'))
end

% brute-force build (this will be costly if most of r is zeros )
%                   (this will be slow if the number of blocks is large)
t = zeros(T,M);
for i=1:M/bc
    for j=1:min(i-1,T)
        t((j-1)*br+1:j*br,(i-1)*bc+[1:bc]) = r(:,(i-j)*bc+1:(i-j+1)*bc);
    end
    t((i-1)*br+1:T,(i-1)*bc+[1:bc]) = c(1:T-br*(i-1),:);
end