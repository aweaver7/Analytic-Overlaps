%%
%analytically calculated top-hat coefficients.
%the waist of basis is w0 and the radius of tophat is rad
%provides fit coefficent (2n,2m) to tophat of unit power, multiplying
%by r*sqrt(pi) gives unit intensity top-hat
function fit_coef = top_hat_fit(n,m,w0,rad)
a0= sqrt(2*factorial(2*n)*factorial(2*m))*w0/rad;
frc = (rad/w0)^2;
ep = exp(-frc);
c0 = 0;
for j=0:n
    for k=0:m
        d0 = ep*sum(sum((frc.^(0:(m+n-(j+k))))./factorial(0:(m+n-(j+k)))));
        d0 = 1-d0;
        f0 = (((-1)^(j+k))/((2^(j+k))*factorial(j)*...
            factorial(k)*factorial(n-j)*factorial(m-k)));
        c0 = c0+ (f0*d0);
    end
end
fit_coef = a0*c0;
end