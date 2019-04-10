%This code presents the analytical solution of the FP
% df/dt = 1/x^2 d/dx (Ca*FunCoef1*df/dx+Cf*FunCoef3*f) 
Q = 2;b = 1;
% ExactFF = @(x,t)(Q*exp(-b*x.^2));
ExactF = @(x,t)(Q*exp(-x.^2));%(x.^2.*(4-x));
source = @(x,t)(x-x); 
ExactQ = @(x,t)(-2*b*x.*Q.*exp(-b.*x.^2));

psi = @(x)(1./x.^2.*(erf(x)-2*x/sqrt(pi).*exp(-x.^2)));

Ca = 1;
Cf = 1;
E = 0;
R = 0;

PDE.term1.Opt = 'Grad';
PDE.term1.FunCoef = @(x)FunCoef1(x);
PDE.term1.Coef =  Cf;

PDE.term2.Opt = 'Diff';
PDE.term2.FunCoef = @(x)FunCoef3(x);
PDE.term2.Coef = Ca;

% E terms
PDE.term3.Opt = 'Grad';
PDE.term3.FunCoef = @(x)(x.^2);
PDE.term3.Coef = -E;

% R terms
PDE.term4.Opt = 'Grad';
PDE.term4.FunCoef = @(x)(x.^3);
PDE.term4.Coef = R;

PDE.BC.q_L = 0; PDE.BC.q_R = 1;
PDE.BC.f_L = 1;  PDE.BC.f_R = 0;

function y = FunCoef1(x)
psi = @(x)(1./x.^2.*(erf(x)-2*x/sqrt(pi).*exp(-x.^2)));
if abs(x)<1e-5
    y = 0;
else
    y = 2*x.^2.*psi(x);
end
end

function y = FunCoef2(x)
psi = @(x)(1./x.^2.*(erf(x)-2*x/sqrt(pi).*exp(-x.^2)));
dpsi = @(x)(2*exp(-x.^2)/sqrt(pi)-(erf(x)-2*x.*exp(-x.^2)/sqrt(pi))./x.^3);
if abs(x)<1e-5
    y = 0;
else
    y = -(1/2)*(erf(x)-2*x.*exp(-x.^2)/sqrt(pi))./x.^2+2*x.*exp(-x.^2)/sqrt(pi);
end
end

function y = FunCoef3(x)
psi = @(x)(1./x.^2.*(erf(x)-2*x/sqrt(pi).*exp(-x.^2)));
if abs(x)<1e-5
    y = 0;
else
    y = x.*psi(x);
end

end