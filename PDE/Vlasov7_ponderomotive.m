function pde=Vlasov7_ponderomotive
% Numerical Example for Vlasov Equation
% This test has given E and non-zero source term

% parameters
k_0=0.5;
A=1;
Lmin=0;Lmax=1;
Vmin=-5;Vmax=+5;
params.TEND = 0.5;

params.k_0 = k_0;
params.A = A;

params.Lmin = Lmin;
params.Lmax = Lmax;
params.Vmin = Vmin;
params.Vmax = Vmax;

pde.Fx_0 = @Fx_0;
pde.Fv_0 = @Fv_0;
pde.Fxv_0 = @Fxv_0;
pde.Ex = @Ex;
pde.Et = @Et;
pde.E = @E;
pde.rho = @rho;
pde.params = params;

pde.solvePoisson = 0;
pde.applySpecifiedE = 1;

pde.checkAnalytic = 0;

% pde.exactE = @exactE;

pde.source1x = @source1x;
pde.source1v = @source1v;
pde.source1t = @source1t;
pde.source1 = @source1;

pde.source2x = @source2x;
pde.source2v = @source2v;
pde.source2t = @source2t;
pde.source2 = @source2;

pde.source3x = @source3x;
pde.source3v = @source3v;
pde.source3t = @source3t;
pde.source3 = @source3;

pde.ExactFx = @ExactFx;
pde.ExactFv = @ExactFv;
pde.ExactFt = @ExactFt;
pde.ExactF = @ExactF;

end

function f=Fx_0(x, params)
% Initial condition for x variable

A = params.A;
k_0 = params.k_0;
Lmax = params.Lmax;
Vmax = params.Vmax;

f=x.*0+1;
end

function f=Fv_0(v, params)
% Initial condition for v variable

A = params.A;
k_0 = params.k_0;
Lmax = params.Lmax;
Vmax = params.Vmax;

f = exp(-params.k_0*v.^2/2);
end
function f=Fxv_0(x,v, params)
A = params.A;
k_0 = params.k_0;
Lmax = params.Lmax;
Vmax = params.Vmax;

f=Fv_0(v).*Fx_0(x);
end
function f=Ex(x,  params)
f=10*exp(-(x-0.5).^2/0.005);;
end
function f=Et(t, params)
f=t*0+1;
end
function f=E(x,t, params)
f=Ex(x).*Et(t);
end
function f=F(x,v,t,  params)
A = params.A;
k_0 = params.k_0;
Lmax = params.Lmax;
Vmax = params.Vmax;

f=t.*x.*(1-x).*(v-Vmax).*(v+Vmax);%.*x.*(Lmax-x);%t*(v-Vmax).*(v+Vmax).*x.*(Lmax-x);
end
function f=rho(x,t, params)
% we do note use rho in this test, so I just write an arbitrary function
A = params.A;
k_0 = params.k_0;
Lmax = params.Lmax;
Vmax = params.Vmax;

f=x-x+1;
end

% function f = exactE(x, params)
% % Exact solution for E
% f=10*exp(-(x-0.5).^2/0.005);
% end
% source term--fully seperable functions
% source = source1+source2+source3
% source term 1
function f = source1x(x,params)
Vmax=5;Lmax=1;
f = x.*0;
end
function f = source1t(t)
Vmax=5;Lmax=1;
f = t.*0;
end
function f = source1v(v,params)
Vmax=5;Lmax=1;
f = v.*0;
end
function f = source1(x,v,t)
Vmax=5;Lmax=1;
f = source1x(x).*source1v(v).*source1t(t);
end

% source term 2
function f = source2x(x,params)
Vmax=5;Lmax=1;
f = x.*0;
end
function f = source2t(t)
Vmax=5;Lmax=1;
f = t.*0;
end
function f = source2v(v,params)
Vmax=5;Lmax=1;
f = v.*0;
end
function f = source2(x,v,t)
Vmax=5;Lmax=1;
f = source2x(x).*source2v(v).*source2t(t);
end

% source term 3
function f = source3x(x,params)
Vmax=5;Lmax=1;
f = x.*0;
end
function f = source3t(t)
Vmax=5;Lmax=1;
f = t.*0;
end
function f = source3v(v,params)
Vmax=5;Lmax=1;
f = v.*0;
end
function f = source3(x,v,t)
Vmax=5;Lmax=1;
f = source3x(x).*source3v(v).*source3t(t);
end

% Exact F
function f=ExactFx(x,params)
Vmax=5;Lmax=1;
f = x.*(1-x);
end
function f=ExactFv(v,params)
Vmax=5;Lmax=1;
f = (v-Vmax).*(v+Vmax);
end
function f=ExactFt(t)
Vmax=5;Lmax=1;
f=t;
end
function f=ExactF(x,v,t)
Vmax=5;Lmax=1;
f = ExactFx(x).*ExactFv(v).*ExactFt(t);
end
