function pde = fokkerplanck1_4p1b
% Problem 4.1b from the RE paper.  
% df/dt == -d/dz ( (1-z^2)f )
%
% Run with
%
% explicit
% fk6d(fokkerplanck1_4p1b,4,2,0.02,[],[],0,[])
%
% implicit
% fk6d(fokkerplanck1_4p1b,6,4,0.5,[],[],1,[],[],1.0)

pde.CFL = 0.01;
sig = 0.1;


%% Setup the dimensions
% 
% Here we setup a 1D problem (x)

    function ret = phi(z,t)
        ret = tanh(atanh(z)-t);
    end
    function ret = f0(z)
        ret = exp(-z.^2/sig^2);
    end
    function ret = soln(z,t)
        p = phi(z,t);
        t1 = 1-p.^2;
        t2 = 1-z.^2;
        t3 = f0(p);
        ret = t1./t2.*t3;
    end

dim_z.BCL = 'D'; % dirichlet
dim_z.BCR = 'D';
dim_z.domainMin = -1;
dim_z.domainMax = +1;
dim_z.init_cond_fn = @(z,p) soln(z,0);

%%
% Add dimensions to the pde object
% Note that the order of the dimensions must be consistent with this across
% the remainder of this PDE.

pde.dimensions = {dim_z};

%% Setup the terms of the PDE
%
% Here we have 1 term1, having only nDims=1 (x) operators.

%% 
% -d/dz ( (1-z^2)*f )

term2_z.type = 'grad'; % grad (see coeff_matrix.m for available types)
term2_z.G = @(z,p,t,dat) -(1-z.^2); % G function for use in coeff_matrix construction.
term2_z.LF = -1; % Upwind 

term2 = {term2_z};

%%
% Add terms to the pde object

pde.terms = {term2};

%% Construct some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

pde.params = params;

%% Add an arbitrary number of sources to the RHS of the PDE
% Each source term must have nDims + 1

%%
% Add sources to the pde data structure
pde.sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

pde.analytic_solutions_1D = { ...
    @(z,p,t) soln(z,t), ...
    @(t,p) 1 
    };

%%
% Function to set time step
    
    function dt=set_dt(pde)
        Lmax = pde.dimensions{1}.domainMax;
        LevX = pde.dimensions{1}.lev;
        CFL = pde.CFL;
        dt = Lmax/2^LevX*CFL;
    end

pde.set_dt = @set_dt;

end

