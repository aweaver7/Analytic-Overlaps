%%
%Function:
%   noll_convert
%Purpose:
%   Convert noll index of zernike polynomial to 2-d zernike indices.
%Syntax:
%   cnv = noll_convert(num)
%Input:
%   num:    Noll zernike index.
%Output:
%   cnv:    Array of converted index [m,n] where m is angular argument and
%           n is the max order.
%
%%
function cnv = noll_convert(num)
n = 0;
while (n+1)*(n+2)/2 < num
    n = n+1;
end 
res = num - (n*(n+1)/2);
pel = abs(mod(n,2)-mod(res,2));
res = res - pel;
b = mod(n,4);
if (pel > 0)
    if (b > 1)
        m = res;
    else
        m = -res;
    end
else
    if (b > 1)
        m = -res;
    else
        m = res;
    end
end
cnv = [m,n];
end