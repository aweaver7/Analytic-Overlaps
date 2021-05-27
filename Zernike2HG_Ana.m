%%
%Function:
%   Zernike2HG_Ana
%Purpose:
%   Analytically calculate fractional amplitude of 
%   beam in desired outgoing HG mode with desired
%   incoming HG mode and Zernike Map term.
%Syntax
%   Fnc_Hndl = Zernike2HG_Ana(w0in,w0out,z0in,z0out,lamIn,lamOut,zApp,rad)
%Inputs:
%   w0in:   Incoming beam waist.
%   w0out:  Outgoing beam waist.
%   z0in:   Incoming beam waist location.
%   z0out:  Outgoing beam waist location.
%   lamIn:  Wavelength of incoming beam.
%   lamOut: Wavelength of outgoing beam.
%   zApp:   Location of the aperture or pupil along the optical axis.
%   rad:    Radius of the aperture or pupil.
%Output:
%   Fnc_Hndl:   Handle to function of 6 variables,
%               Fnc_Hndl(a,b,c,d,e,f) gives amount
%               of incoming mode (e,f) leaving in
%               mode (a,b) after going through zernike
%               (c,d) where c is actually the upper
%               (smaller or angular) index
%Changes Tracked:
%   changed r dependence from kappa to belle, 1/r^d was in kappa r^2c in
%   belle, changed to r^(2C-d) in belle.
%   fixed guoy phase term to include initial offset
%%
function Fnc_Hndl = Zernike2HG_Ana(w0in,w0out,z0in,z0out,lamIn,lamOut,zApp,rad)
%1 indices count for outgoing, 2 for incoming
w2 = w0in;
w1 = w0out;
dz2 = zApp - z0in;
dz1 = zApp - z0out;
k2 = 2*pi()/lamIn;
k1 = 2*pi()/lamOut;
zr2 = pi()*(w2^2)/lamIn;
zr1 = pi()*(w1^2)/lamOut;
fr1 = dz1./zr1;
fr2 = dz2./zr2;
r = rad;
W1 = @(a) ((w1.^2)+((dz1.*lamOut./(pi().*w1)).^2)).^(a/2);
W2 = @(a) ((w2.^2)+((dz2.*lamIn./(pi().*w2)).^2)).^(a/2);
gam = ((1-(1i.*fr1))./W1(2))+((1+(1i.*fr2))./W2(2));
gu1 = atan(fr1);
gu2 = atan(fr2);
exga = exp(-gam.*(r.^2));
fac = @(a,b,c,d) factorial(a).*factorial(b).*factorial(c).*factorial(d);
kap = @(a,b,c,d,f,g) factorial(c).*sqrt(fac(a,b,f,g)).*(2.^(a+b+f+g+1)).*exp(1i.*((k1.*dz1)-(k2.*dz2)+...
    ((f+g+1).*gu2)-((a+b+1).*gu1)))./(pi().*W1(a+b+1).*W2(f+g+1).*(gam.^((a+b+d+f+g+2)/2)));
bel = @(a,b,c,d,f,g,A,B,C,D,F,G,K) ((-gam./8).^(A+B+F+G)).*((-1).^(C+D)).*((gam).^(C)).*...
    W1(2.*(A+B)).*W2(2.*(F+G)).*factorial(((a+b+d+f+g)./2)-(A+B+C+F+G)).*...
    (r.^((2.*C)-d)).*factorial(d-C).*gamma(((b+g+1+K)./2)+D-(B+G)).*gamma(((a+c+f+1-K)./2)-(A+D+F))./...
    (factorial(C).*fac(A,B,F,G).*factorial(((a+b+c+f+g)/2)-(A+B+F+G)).*...
    fac(a-(2.*A),b-(2.*B),f-(2.*F),g-(2.*G)).*fac((((d+c)./2)-C),(((d-c)./2)-C),((2.*D)+K),(c-(2.*D)-K)));
if isinf(r)
    exPa = @(a,b,c,d,f,g,A,B,C,D,F,G) 1;
else
    exPa = @(a,b,c,d,f,g,A,B,C,D,F,G) 1-(exp(-gam.*(r.^2)).*...
        sum(((gam.*(rad.^2)).^(0:(((a+b+d+f+g)./2)-(A+B+C+F+G))))./factorial(0:(((a+b+d+f+g)/2)-(A+B+C+F+G)))));
end

Fnc_Hndl = @(a,b,c,d,f,g) ovlpfnc(a,b,c,d,f,g);
    function Val = ovlpfnc(a,b,c,d,f,g)
        ma = numel(a);
        mb = numel(b);
        mc = numel(c);
        md = numel(d);
        mf = numel(f);
        mg = numel(g);
        Val = zeros(ma,mb,mc,md,mf,mg);
        for Ma = 1:ma
            aa = a(Ma);
            for Mb = 1:mb
                bb = b(Mb);
                for Mc = 1:mc
                    Tc = mod(c(Mc),2);
                    cc = c(Mc);
                    for Md = 1:md
                        Td = mod(d(Md),2);
                        dd = d(Md);
                        for Mf = 1:mf
                            ff = f(Mf);
                            for Mg = 1:mg 
                                gg = g(Mg);
                                if (abs(cc) > aa+bb+ff+gg)
                                    Val(Ma,Mb,Mc,Md,Mf,Mg) = 0;
                                    tk = 0;
                                elseif abs(cc) > abs(dd)
                                    Val(Ma,Mb,Mc,Md,Mf,Mg) = 0;
                                    tk = 0;
                                elseif Tc ~= Td
                                    Val(Ma,Mb,Mc,Md,Mf,Mg) = 0;
                                    tk = 0;
                                elseif (cc < 0) && ...
                                        ((mod(bb+gg,2)~=1) ||...
                                        (mod(aa+ff+dd,2)~=1))
                                    Val(Ma,Mb,Mc,Md,Mf,Mg) = 0;
                                    tk = 0;
                                elseif (cc >=0) && ((mod(bb+gg,2)~=0) ||...
                                        (mod(aa+ff+dd,2)~=0))
                                    Val(Ma,Mb,Mc,Md,Mf,Mg) = 0;
                                    tk = 0;
                                elseif (cc < 0)
                                    tl = 1;
                                    tk = 1;
                                elseif (cc >= 0)
                                    tl = 0;
                                    tk = 1;
                                end
                                if tk == 1
                                    ccc = abs(cc);
                                    for A = 0:fix(aa/2)
                                        for B = 0:fix(bb/2)
                                            for C = 0:fix((dd-ccc)/2)
                                                for D = 0:fix((ccc-tl)/2)
                                                    for F = 0:fix(ff/2)
                                                        for G = 0:fix(gg/2)
                                                            Val(Ma,Mb,Mc,Md,Mf,Mg) = ...
                                                                Val(Ma,Mb,Mc,Md,Mf,Mg) + ...
                                                                (kap(aa,bb,ccc,dd,ff,gg)*exPa(aa,bb,ccc,dd,ff,gg,A,B,C,D,F,G)*bel(aa,bb,ccc,dd,ff,gg,A,B,C,D,F,G,tl));
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end