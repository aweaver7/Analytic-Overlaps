%%
%Function: TEM_Analytic_Overlap
%Purpose:   Return function giving overlap of TEM modes on
%           an apperture
%Syntax:
%   over_fnc = TEM_overlap(w1,w2,z1,z2,l,z,r)
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
%               having waist w2 at z2 that leaves in after being clipped
%               by an apperture of radius r at position z in the mode
%               j,k of waist size w1 at z1.
%           
%Inputs:
%   w1: The waist of the output beam.
%   w2: The waist of the input beam.
%   z1: The waist position of the output beam.
%   z2: The waist position of the input beam.
%   l: The wavelength of the beams.
%   z: The z-coordinate of the apperture.
%   r:  The radius of the apperture.
%%%
function over_fnc=TEM_Analytic_Overlap(w1,w2,z1,z2,l,z,r)
raj = w1;
rbj = w2;
w1 = @(a) raj.^a;
w2 = @(a) rbj.^a;
R = @(a) r^a;
k = 2*pi()/l;
zr1 = pi()*w1(2)/l;
zr2 = pi()*w2(2)/l;
dz1 = z-z1;
dz2 = z-z2;
P1P = @(a) (1+(1i*dz1/zr1)).^a;
P2P = @(a) (1+(1i*dz2/zr2)).^a;
P1M = @(a) (1-(1i*dz1/zr1)).^a;
P2M = @(a) (1-(1i*dz2/zr2)).^a;
sig = @(a) ((P1P(-1)*w1(-2))+(P2M(-1)*w2(-2))).^a;
wsr = @(a) w1(a).*w2(a);
W1 = @(a) w1(a).*P1P(a/2).*P1M(a/2);
W2 = @(a) w2(a).*P2P(a/2).*P2M(a/2);
WSR = @(a) W1(a).*W2(a);
delK = k*(z2-z1);
exf = @(a) 1-(exp(-sig(1)*R(2))*sum(((sig(1)*R(2)).^(0:a))./factorial(0:a)));
coef1 = @(a,b,c) exp(1i*delK).*wsr(c+1).*P1M(a+c+1).*P2P(b+c+1).*...
    P1P(-a).*P2M(-b)./(2^(a+b+c-1));
coef2 = @(a,b,c,d) sqrt(factorial(a).*factorial(b).*...
    factorial(c).*factorial(d));
coef3 = @(c1,c2,d1,d2,m1,m2,j1,j2,l1,l2) W1(2.*(c2+d2)).*W2(2.*(c1+d1)).*...
    ((-1).^(m1+m2+j1+j2-(c1+c2+d1+d2))).*prod(1:2:(2.*(d1+d2+l2))).*...
    prod(1:2:(2.*(c1+c2+l1))).*(2^(2.*(c1+c2+d1+d2+l1+l2)))./...
    (factorial(m1-c1).*factorial(m2-c2).*factorial(j1-d1).*...
    factorial(j2-d2).*factorial((2.*c1)+l1).*...
    factorial((2.*c2)+l1).*factorial((2.*d1)+l2).*...
    factorial((2.*d2)+l2));
aexf = @(a) sig(-(a+1)).*WSR(-2.*(a+1)).*exf(a);
%up until here we've been setting up coefficients that come up in the 
%analytic expression.
over_fnc = @(a,b,c,d) retove(a,b,c,d);

%this function builds the analytic expression using the terms prepared
%above and returns the function handle we desire with the overlap.
function ovr = retove(i1,i2,i3,i4)
m1 = fix(i1/2);
l1 = fix(i1)-(2*m1);
m2 = fix(i3/2);
h1 = fix(i3)-(2*m2);
j1 = fix(i2/2);
l2 = fix(i2)-(2*j1);
j2 = fix(i4/2);
h2 = fix(i4)-(2*j2);
n1 = numel(i1);
n2 = numel(i2);
n3 = numel(i3);
n4 = numel(i4);
ovr = zeros(n1,n2,n3,n4);
for A = 1:n1
    M1 = m1(A);
    L1 = l1(A);
    for B = 1:n2
        J1 = j1(B);
        L2 = l2(B);
        S1 = M1+J1;
        for C = 1:n3
            M2 = m2(C);
            H1 = h1(C);
            for D=1:n4
                J2 = j2(D);
                H2 = h2(D);
                if (L1 == H1)&&(L2 == H2)
                    mum = 0;
                    S2 = M2+J2;
                    L = L1+L2;
                    for c1=0:M1
                        for c2= 0:M2
                            for d1=0:J1
                                for d2 = 0:J2
                                    mum=mum+...
                                        (coef3(c1,c2,d1,d2,...
                                        M1,M2,J1,J2,L1,L2).*...
                                        aexf(c1+c2+d1+d2+L));
                                end
                            end
                        end
                    end
                    ovr(A,B,C,D) = mum.*coef1(S1,S2,L).*...
                        coef2(i1(A),i2(B),i3(C),i4(D));
                end
            end
        end
    end
end
end
end
 