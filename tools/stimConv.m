function y = stimConv(X,F)
[n,nk] = size(X); 
[ntau,nkf]=size(F);
if nkf~=nk
    error('X and F must have same number of columns');
end

y = zeros(n+ntau-1,1);
for i = 1 : nk 
    xfconv = conv(X(:,i),F(:,i));
    y = y + xfconv;
end
y = y(1:n);
end

