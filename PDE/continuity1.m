function pde = continuity1
% 1D test case using continuity equation, i.e., 
%
% df/dt == -df/dx
%
% Run with
%
% explicit
% asgard(continuity1)
% asgard(continuity1,'lev',4,'deg',3)
%
% implicit
% asgard(continuity1,'timestep_method','CN')
% asgard(continuity1,'timestep_method','CN','CFL',0.1)

%% Setup the dimensions
% 
% Here we setup a 1D problem (x)

dim_x.domainMin = -1;
dim_x.domainMax = +1;
dim_x.init_cond_fn = @(x,p,t) x.*0;

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_x};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE
%
% Here we have 1 term1, having only nDims=1 (x) operators.

%% 
% -df/dx

g1 = @(x,p,t,dat) x.*0-1;
pterm1 = GRAD(num_dims,g1,0,'P','P');

term1_x = TERM_1D({pterm1});
term1   = TERM_ND(num_dims,{term1_x});

%%
% Add terms to the pde object

pde.terms = {term1};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%% Add an arbitrary number of sources to the RHS of the PDE
% Each source term must have nDims + 1 (which here is 2+1 / (x,v) + time) functions describing the
% variation of each source term with each dimension and time.
% Here we define 3 source terms.

%%
% Source 1
s1x = @(x,p,t) cos(2*pi*x);
s1t = @(t,p) cos(t);
source1 = {s1x,s1t};

%%
% Source 2
s2x = @(x,p,t) sin(2*pi*x);
s2t = @(t,p) -2*pi*sin(t);
source2 = {s2x,s2t};

%%
% Add sources to the pde data structure
pde.sources = {source1,source2};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

a_x = @(x,p,t) cos(2*pi*x);
a_t = @(t,p) sin(t);

pde.analytic_solutions_1D = {a_x,a_t};

%%
% Function to set time step

    function dt=set_dt(pde,CFL)
        
        dim = pde.dimensions{1};
        lev = dim.lev;
        xMax = dim.domainMax;
        xMin = dim.domainMin;
        xRange = xMax-xMin;
        dx = xRange/(2^lev);
        dt = CFL*dx;
    end

pde.set_dt = @set_dt;

end


