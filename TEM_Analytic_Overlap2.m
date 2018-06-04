%%
%Function: TEM_Analytic_Overlap2
%Purpose:   Return function giving overlap of TEM modes on
%           an apperture
%Syntax:
%   over_fnc = TEM_Analytic_Overlap2(w1,w2,z1,z2,l1,l2,z,r)
%Output:    
%   over_fnc:   The overlap function of HG-TEM mode with  
%               the conjugate of a HG-TEM mode, should
%               be be a 4 argument function so that
%               over_fnc(a,b,c,d) gives the overlap
%               of mode (c,d) with the conjugate of
%               (a,b) over the desired sized apperture
%               with first (conjugate) beam having desired
%               properties (w1,z1,l1) and the second
%               having corresponding properties for real
%               (incoming) beam. In other words the function gives
%               the amount of a beam with waist w2 at location
%               z2 clipped by an apperture of radius r at position
%               z leaving in a mode with waist w1 located at z1.
%               over_fnc(j,k,l,m) gives the amount of a TEM l,m mode
%               having waist w2 at z2 with wavelength l2 that leaves after 
%               being clipped by an apperture of radius r at position z in 
%               the mode j,k of waist size w1 at z1 with wavelength l1.
%           
%Inputs:
%   w1: The waist of the output beam.
%   w2: The waist of the input beam.
%   z1: The waist position of the output beam.
%   z2: The waist position of the input beam.
%   l1: The wavelength of the outgoing beam.
%   l2: The wavelength of the incoming beam.
%   z: The z-coordinate of the apperture.
%   r:  The radius of the apperture.
%%%
function over_fnc=TEM_Analytic_Overlap2(w1,w2,z1,z2,l1,l2,z,r)
%setting up wavenumber of incoming and outgoing beam
ki = 2*pi()/l2;
ko = 2*pi()/l1;
%defining distance of apperture from incoming and outgoing beams
dzi = z-z2;
dzo = z-z1;
%setting up rayleigh range of incoming and outgoing
zri = ki*(w2^2)/2;
zro = ko*(w1^2)/2;
Mi = @(a) (w2.*(1-(1i*dzi/zri))).^a;
Mo = @(a) (w1.*(1-(1i*dzo/zro))).^a;
Pi = @(a) (w2.*(1+(1i*dzi/zri))).^a;
Po = @(a) (w1.*(1+(1i*dzo/zro))).^a;
%phase along optical path
exp1 = (ko.*dzo)-(ki.*dzi);
%relabeling wo and wi for easier use
wo = w1;
wi = w2;
%for the second exponent term
exp2 = @(a) ((r.^2).*(((wi*Mi(1)).^(-1))+((wo*Po(1)).^(-1)))).^a;
%now creating factors for terms raised to f+g and h+l
fac1 = @(a) ((-Mo(1)/(2*wo))*(1+(wo*Po(1)/(wi*Mi(1))))).^a;
fac2 = @(a) ((-Pi(1)/(2*wi))*(1+(wi*Mi(1)/(wo*Po(1))))).^a;
%building actual exponents
EX1 = exp(i*exp1);
EX2 = exp(-exp2(1));
%now creating factor for in front of summations besides
%first exponent.
%first for total factorial
totfac = @(a,b,c,d) factorial(a).*factorial(b).*...
    factorial(c).*factorial(d);
%note here (a,b,c,d) corresponds to (m,n,c,d)
bfac1 = @(a,b,c,d) 2*sqrt(totfac(a,b,c,d)./...
    (((Po(1)/wi)+(Mi(1)/wo)).^(a+b+c+d+2))).*...
    sqrt(Po(c+d-(a+b)).*Mi(a+b-(c+d))).*EX1;
%here is term for first set of summations, (f,g,h,l)
%corresponds to same terms, doesn't include factorials involving
%initial indices
afac1 = @(f,g,h,l) fac1(f+g).*fac2(h+l)./totfac(f,g,h,l);
%now for total term that is argument of final index
afac2 = @(k) (EX2./factorial(k)).*exp2(k);
over_fnc = @(m,n,c,d) tmp_fnc_tot(m,n,c,d);
function outp = tmp_fnc_tot(A,B,C,D)
n1 = numel(A);
n2 = numel(B);
n3 = numel(C);
n4 = numel(D);
indm = floor(A/2);
indn = floor(B/2);
indc = floor(C/2);
indd = floor(D/2);
Amod2 = mod(A,2);
Bmod2 = mod(B,2);
Cmod2 = mod(C,2);
Dmod2 = mod(D,2);
outp = zeros(n1,n2,n3,n4);
for ind1 = 1:n1
    for ind2 = 1:n2
        for ind3 = 1:n3
            for ind4 = 1:n4
                if (Amod2(ind1)~=Cmod2(ind3))||...
                        (Bmod2(ind2)~=Dmod2(ind4))
                    outp(ind1,ind2,ind3,ind4) = 0;
                else
                    m = A(ind1);
                    mind = indm(ind1);
                    n = B(ind2);
                    nind = indn(ind2);
                    c = C(ind3);
                    cind = indc(ind3);
                    d = D(ind4);
                    dind = indd(ind4);
                    TT = 0;
                    for f = 0:mind
                        for g = 0:nind
                            for h = 0:cind
                                for l = 0:dind
                                    FC = sum(sum(afac2(0:(((m+n+c+d)/2)-(f+g+h+l)))));
                                    FCC = 1-FC;
                                    FCCC = afac1(f,g,h,l).*...
                                        (factorial(n+d-(2*(g+l)))./...
                                        (factorial(n-(2*g)).*factorial(d-(2*l)))).*...
                                        (factorial(m+c-(2*(f+h)))./...
                                        (factorial(m-(2*f)).*factorial(c-(2*h))))./...
                                        (factorial(((n+d)/2)-(g+l)).*...
                                        factorial(((m+c)/2)-(f+h))).*FCC;
                                    TT = TT+FCCC;
                                end
                            end
                        end
                    end
                    TT = TT.*bfac1(m,n,c,d);
                    outp(ind1,ind2,ind3,ind4) = TT;
                end
            end
        end
    end
end
end
end