%%
%Program:
%   Field_Act
%Purpose:
%   Render the field over some region.
%Syntax:
%   Fnc_Hndl = Field_Act(coef_vec,waist,waist_location,lambda,zoff)
%Inputs:
%   coef_vec:       Coefficient vector of HG mode coefficients in simtools
%                   format.
%   waist:          Waist size
%   waist_location: Location of the waist.
%   lambda:         Wavelength of the beam.
%   zoff:           Offset of z to use for phase in k*z part of phase,
%                   uses k*(zoff) instead of k*z, e^(-ik(zoff))
%Outputs:
%   Fnc_Hndl:   Handle to function of 3 variables (x,y,z) that gives the
%               phase at these points, x,y, and z must either all be the
%               same size or singletons.
%%
function Fnc_Hndl = Field_Act(coef_vec,waist,waist_location,lambda,zoff)
arr = coef_vec;
w0 = waist;
z0 = waist_location;
lam = lambda;
zr = pi()*w0*w0/lam;
k = 2*pi()/lam;
psi = @(z) atan2((z-z0),zr);
w2 =  @(z) ((w0*w0)+((z-z0).*(z-z0).*lam.*lam./pi()./pi()./w0./w0));
%hv gives value of m+n, argument is index
hv = @(m) ceil((sqrt(9+(8*(m-1)))-3)/2);
%nv gives value of n, argument is index
nv = @(m) hv(m)+m-((hv(m)+1).*(hv(m)+2)/2);
bls = 0.*arr;
bls(:) = 1:numel(arr);
h = hv(bls);
n = nv(bls);
m = h-n;
q = @(z) 1+(1i.*(z-z0)./zr);
Q = @(z) conj(q(z));
cnm = arr./sqrt((2.^h).*factorial(m).*factorial(n));
rh = @(x,y) (x.^2)+(y.^2);
Fnc_Hndl = @(x,y,z) exp(-1i.*k.*(zoff)).*(sqrt(2/pi())./(w0.*Q(z))).*exp(-(rh(x,y)./((w0^2).*Q(z)))).*smtot(sqrt(2./w2(z)).*x,sqrt(2./w2(z)).*y,q(z),Q(z),cnm,m,n,h);
end
%now for the smtot function
function val = smtot(a,b,q,Q,cnm,m,n,h)
if (numel(a) > 1)
    b = (a.*0) + b;
    q = (a.*0) + q;
    Q = (a.*0) + Q;
elseif (numel(b) > 1)
    a = (b.*0)+a;
    q = (b.*0)+q;
    Q = (b.*0)+Q;
else
    a = (q.*0)+a;
    b = (q.*0)+b;
end
val = a.*0;
Hx = @(t) polyval(hrmcrt(t),a);
Hy = @(t) polyval(hrmcrt(t),b);
for mm = 1:numel(cnm)
    val = val + (cnm(mm).*Hx(m(mm)).*Hy(n(mm)).*((q./Q).^(h(mm)./2)));
end
end


%hermite function creator
function hk = hrmcrt(n)
if n==0 
    hk = 1;
elseif n==1
    hk = [2 0];
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