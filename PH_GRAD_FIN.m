%%
%Function: 
%   PH_GRAD_FIN
%Purpose:
%   Function handle to analytic calculation of phase
%   derivative of desired beam over desired region.
%Syntax:
%   [Fnc_Hnd, suc] = PH_GRAD_FIN(wst,wst_loc,lambda,...
%                   coef_vec,xreg,yreg,zreg)
%Inputs:
%   wst:        Waist size of basis to build beam from.
%   wst_loc:    Location of waist in basis used.
%   lambda:     Wavelength of basis.
%   coef_vec:   Vector of coefficients for HG mode in same
%               format as would be used for simtools.
%   xreg:       X-region to evaluate over.
%   yreg:       Y-region to evaluate over.
%   zreg:       Z-region to evaluate over.
%Outputs:
%   Fnc_Hnd:    Handle to function of 3 values describing
%               decomposition of derivative into gradient
%               components, Fnc_Hnd(a,b,c) returns the 
%               derivative along ae_x + be_y + ce_z.
%   suc:        1 or 0 depending on whether it was
%               successful or not.
%%
function [Fnc_Hnd, suc] = PH_GRAD_FIN(wst,wst_loc,lambda,coef_vec,xreg,yreg,zreg)
suc = 0;
w0 = wst;
z0 = wst_loc;
z = zreg-z0;
k = 2*pi()/lambda;
zr = pi()*(w0^2)/lambda;
x = xreg;
y = yreg;
nz = numel(z);
nx = numel(x);
ny = numel(y);
%setting sizes same or trying of field
if nx > 1
    y = y.*ones(size(x));
    z = z.*ones(size(x));
elseif ny > 1
    x = x.*ones(size(y));
    z = z.*ones(size(y));
elseif nz > 1
    x = x.*ones(size(z));
    y = y.*ones(size(z));
end
Sdim = size(x);
rh2 = (x.^2)+(y.^2);
frc = z./zr;
w2 = (w0^2).*(1+(frc.^2));
xder = -2.*x.*frc./w2;
yder = -2.*y.*frc./w2;
zderwoK = (-rh2./zr./w2)+(2*(w0^2).*rh2.*(frc.^2)./zr./w2./w2)+((w0.^2)./zr./w2);
Sdim = size(xder);
x = sqrt(2).*x./w0./sqrt(1+(frc.^2));
y = sqrt(2).*y./w0./sqrt(1+(frc.^2));
ps = atan2(z,zr);
armg = abs(coef_vec);
arph = angle(coef_vec);
%now for rendering other parts of derivative.
for M = 1:Sdim(1)
    %these terms come from the latex variables of same name
   [T,D,X1,X2,Y1,Y2,P1,P2] = DERV(armg,arph,x(M,:),y(M,:),w2(M,:),ps(M,:));
   xyoutfc = 2.*sqrt(2./w2(M,:));
   xder(M,:) = xder(M,:) + (xyoutfc.*((D.*X1)-(T.*X2))./((T.^2)+(D.^2)));
   yder(M,:) = yder(M,:) + (xyoutfc.*((D.*Y1)-(T.*Y2))./((T.^2)+(D.^2)));
   outz1 = 2.*z(M,:).*w0.*w0./zr./zr./w2(M,:);
   outz2 = w0.*w0./zr./w2(M,:);
   zderwoK(M,:) = zderwoK(M,:) + (outz1.*((x(M,:).*((D.*X1)-(T.*X2)))+(y(M,:).*((D.*Y1)-(T.*Y2))))./((T.^2)+(D.^2)))+(outz2.*((D.*P1)+(T.*P2))./((T.^2)+(D.^2)));
end
Fnc_Hnd = @(a,b,c) ((a.*xder)+(b.*yder)+(c.*zderwoK)) - (c.*k);
suc = 1;
end

function [T,D,X1,X2,Y1,Y2,P1,P2] = DERV(armg,arph,x,y,w2,ps)
%hv gives value of m+n, argument is index
hv = @(m) ceil((sqrt(9+(8*(m-1)))-3)/2);
%nv gives value of n, argument is index
nv = @(m) hv(m)+m-((hv(m)+1)*(hv(m)+2)/2);
%mv gives value of m, argument is index
mv = @(m) hv(m)-nv(m);
%First zeroing each of the terms
Hx = @(n) polyval(hrmcrt(n),x);
Hy = @(n) polyval(hrmcrt(n),y);
T = 0.*x;
D = T;
X1= T;
X2= T;
Y1= T;
Y2= T;
P1= T;
P2= T;
%creating mag function
cv = @(m) armg(m)./sqrt((2.^hv(m)).*factorial(mv(m)).*factorial(nv(m)));
%creating cosine or sine term
snv = @(m) sin((hv(m).*ps)+arph(m));
csv = @(m) cos((hv(m).*ps)+arph(m));
%combining these for space
tsv = @(m) cv(m).*snv(m);
tcv = @(m) cv(m).*csv(m);
numm = numel(armg);
for m = 1:numm
    T = T + (tsv(m).*Hx(mv(m)).*Hy(nv(m)));
    D = D + (tcv(m).*Hx(mv(m)).*Hy(nv(m)));
    X1 = X1 + (mv(m).*tsv(m).*Hx(mv(m)-1).*Hy(nv(m)));
    X2 = X2 + (mv(m).*tcv(m).*Hx(mv(m)-1).*Hy(nv(m)));
    Y1 = Y1 + (nv(m).*tsv(m).*Hx(mv(m)).*Hy(nv(m)-1));
    Y2 = Y2 + (nv(m).*tcv(m).*Hx(mv(m)).*Hy(nv(m)-1));
    P1 = P1 + (hv(m).*tcv(m).*Hx(mv(m)).*Hy(nv(m)));
    P2 = P2 + (hv(m).*tsv(m).*Hx(mv(m)).*Hy(nv(m)));
end
end



function hk = hrmcrt(n)
if n==0 
    hk = 1;
elseif n==1
    hk = [2 0];
elseif n<0
    hk = 0;
else    
    hkm2 = zeros(1,n+1);
    hkm2(n+1) = 1;
    hkm1 = zeros(1,n+1);
    hkm1(n) = 2;
    for k=2:n        
        hk = zeros(1,n+1);
        for e=n-k+1:2:n
            hk(e) = 2*(hkm1(e+1) - (k-1)*hkm2(e));
        end       
        hk(n+1) = -2*(k-1)*hkm2(n+1);        
        if k<n
            hkm2 = hkm1;
            hkm1 = hk;
        end        
    end    
end
end