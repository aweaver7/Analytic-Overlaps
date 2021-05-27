%%
%Function:
%   Zernike_Phase
%Purpose:
%   To create a phase map given a set of coefficients in Noll
%   sequential indexing.
%Syntax:
%   Version 1:
%       [Phase_Fnc, success] = Zernike_Phase(r,coef_vector)
%   Version 2:
%       [Phase_Fnc, success] = Zernike_Phase(r,c1,c2,...)
%Inputs:
%   coef_vector: 1-d vector of zernike coefficients.
%   c1-cN:  Zernike coefficients written sequentially.
%Output:
%   Phase_Fnc:  Handle to function for phase over the disc. 2-d
%               function of (rho,phi)
%   success:    0 if it was unsuccessful for some reason, 1 if it worked.
%%
function [Phase_Fnc, success] = Zernike_Phase(r,varargin)
success = 0;
Phase_Fnc = @(rho,phi) 0;
if nargin > 2
     nums = [];
     for l = 1:(nargin-1)
         if numel(varargin{l}) > 1
             return
         else
             nums = [nums,varargin{l}];
         end
     end
else
    nums = cell2mat(varargin);
end
nums = nums(:);
lng = numel(nums);
success = 1;
for l = 1:lng
    arr = noll_convert(l);
    bll = Zernike_Create(arr(1),arr(2),r);
    Phase_Fnc = @(rho,phi) Phase_Fnc(rho,phi) + (nums(l)*bll(rho,phi));  
end
end

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