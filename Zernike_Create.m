%%
%Function: Zernike_Create
%Purpose:
%   Creates a handle to the zernike polynomial of desired order.
%Syntax:
%   [fnc_hndl,suc] = Zernike_Create(index1,index2,radius)
%Inputs:
%   index1: The first (or upper) index of the zernike polynomial.
%   index2: The second (or lower) index of the zernike polynomial.
%   radius: The radius of the disc over which the polynomials
%           are to be orthonormal.
%Outputs:
%   fnc_hndl:   Handle to the zernike polynomial, it is a function of
%               (rho,theta) for the 2 variables.
%   suc:        Whether it was successful.
%%
function [fnc_hndl,suc] = Zernike_Create(index1,index2,radius)
m = index1;
n = index2;
r = radius;
if m ~= fix(m)
    fnc_hndl = @(rho,theta) 0;
    suc = 0;
    return
elseif n ~= fix(n)
    fnc_hndl = @(rho,theta) 0;
    suc = 0;
    return
elseif m < 0
    M = abs(m);
    fncp1 = @(theta) sin(M*theta);
    if n < M
        fnc_hndl = @(rho,theta) 0;
        suc = 0;
        return
    end
else
    M = m;
    fncp1 = @(theta) cos(M*theta);
    if n < M
        fnc_hndl = @(rho,theta) 0;
        suc = 0;
        return
    end
end
fncp2 = @(rho) 0;
if mod(n,2) ~= mod(M,2)
    fncp2 = @(rho) 0;
else
    binc = @(a,b) factorial(a)/(factorial(b)*factorial(a-b));
    for k = 0:((n-M)/2)
        fncp2 = @(rho) fncp2(rho) + (((-1)^k)*binc(n-k,k)*...
            binc(n-(2*k),(((n-m)/2)-k))*((rho/r).^(n-(2*k))));
    end
end
suc = 1;
fnc_hndl = @(rho,phi) fncp1(phi).*fncp2(rho);
end