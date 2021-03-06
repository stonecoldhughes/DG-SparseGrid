function pde = projecti_diff1
% Example PDE using the 1D Diffusion Equation. This example PDE is
% time dependent (although not all the terms are time dependent). This
% implies the need for an initial condition. 
% PDE:
% 
% df/dt = nu * d^2 f/dx^2
%
% Domain is [0,2]
%
%
% d^2 f / dx^2 becomes
%
% dq/dx with free (neumann BC)
%
% and
%
% q=df/dx with Dirichlet BCs specified by analytic solution
%
% Run with
%
% explicit
% asgard(projecti_diff1);
%
% implicit
% asgard(projecti_diff1,'timestep_method','CN');
%
% FE
% asgard(projecti_diff1,'lev',4,'deg',4,'timestep_method','FE', 'dt',0.01,'num_steps',16)

pde.CFL = 0.01;

%% Setup the dimensions
% 
% Here we setup a 1D problem (x,y)

k = pi/2;
nu = 0.01;
soln_x = @(x) sin(k*x);
soln_t = @(t) exp(-nu*k^2*t);

dim_x.name = 'x';
dim_x.domainMin = 0;
dim_x.domainMax = 2;
dim_x.init_cond_fn = @(x,p,t) soln_x(x)*soln_t(t);


%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_x};
num_dims = numel(pde.dimensions);

%% Setup the terms of the PDE
%
% Here we have 1 term, with each term having nDims (x and y) operators.

%% 
% Setup the d^2_dx^2 term

% term1
%
% eq1 :  df/dt   == d/dx g1(x) q(x,y)   [grad,g1(x)=1, BCL=N, BCR=D]
% eq2 :   q(x,y) == d/dx g2(x) f(x,y)   [grad,g2(x)=1, BCL=D, BCR=D]
%
% coeff_mat = mat1 * mat2

g1 = @(x,p,t,dat) x.*0+nu;
g2 = @(x,p,t,dat) x.*0+1;

pterm1 = GRAD(num_dims,g1,+1,'N','N');
pterm2 = GRAD(num_dims,g2,-1,'D','D');

term1_x = TERM_1D({pterm1,pterm2});
term1   = TERM_ND(num_dims,{term1_x});

%%
% Add terms to the pde object

 pde.terms = {term1};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%%
% Sources

pde.sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...
    @(x,p,t) soln_x(x), ...
    @(t,p) soln_t(t) 
    };

    function dt=set_dt(pde,CFL)
        
        dims = pde.dimensions;
        
        % for Diffusion equation: dt = C * dx^2
        
        lev = dims{1}.lev;
        dx = 1/2^lev;
        dt = CFL*dx^2;
        
    end

pde.set_dt = @set_dt;

end

%%
% Function to set time step
