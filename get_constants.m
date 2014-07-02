function [cst1,cst2,cst3,cst4,cst5,cst6,cst7,cst8,...
          cst9,cst10,cst11] = get_constants(beta,gamma,c)

    
%c = 299792.458;  % speed of light in km/s 

cst1  = 2*(beta+gamma)/c^2;
cst2  = (2*beta-1)/c^2;
cst3  = gamma/c^2;
cst4  = (1+gamma)/c^2;
cst5  = 2*(1+gamma)/c^2;
cst6  = 3/(2*c^2);
cst7  = 1/(2*c^2);
cst8  = 1/c^2;
cst9  = double(2+2*gamma);
cst10 = double(1+2*gamma);
cst11 = (3+4*gamma)/(2*c^2);
      

end