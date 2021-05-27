%%
%Program:
%   Coef_Frmt_Zern
%Purpose:
%   Put coefficients in an array that SIMTOOLS can use.
%Syntax:
%   coef_array = Coef_Frmt_Zern(w0,w1,z0,z1,lam0,lam1,z,rad,ind1,ind2,num,arr);
%Output:
%   coef_array: Array of Coefficients in order that is
%               required by the SIMTOOLS
%Inputs:
%   w0:     Waist size of incoming basis.
%   w1:     Waist size of outgoing beam.
%   z0:     Waist location of incoming basis.
%   z1:     Waist location of outgoing beam.
%   lam0:   Wavelength of incoming beam.
%   lam1:   Wavelength of outgoing beam.
%   z:      Location of Aperture
%   rad:    Radius of Aperture.
%   ind1:   First index of incoming mode.
%   ind2:   Second index of incoming mode.
%   num:    Max number indices add up to.
%   arr:    Array of Zernike Coefficients for phase map with beam
%%
function coef_array = Coef_Frmt_Zern(w0,w1,z0,z1,lam0,lam1,z,rad,ind1,ind2,num,arr)
coffnc = Zernike2HG_Ana(w0,w1,z0,z1,lam0,lam1,z,rad);
coffnc = @(m1,m2,m3,m4) coffnc(m1,m2,m3,m4,ind1,ind2);
coef_array = [];
for n1 = 0:num
    for n2 = 0:n1
        coef_array = [coef_array,coffnc(n1-n2,n2,0,0)];
    end
end
for n1 = 1:numel(arr)
    ticc = 1;
    for n2 = 0:num
        for n3 = 0:n2
            nds = noll_convert(n1);
            m = nds(1);
            n = nds(2);
            coef_array(ticc) = coef_array(ticc)+(1i.*arr(n1).*coffnc(n2-n3,n3,m,n));
            ticc = ticc+1;
        end
    end
end
end
