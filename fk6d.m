function [fval,err] = fk6d(pde,Lev,Deg,TEND,quiet,compression,implicit)

%% MATLAB (reference) version of the DG-SG solver
% The default execution solves the Vlasov-Poisson system of equations
%
% $$f_t + v\frac{\partial f}{\partial x} + E\left(x,t\right)\frac{\partial
% f}{\partial v} f=0$$
%
% $$-\frac{\partial^2\phi}{\partial x^2} = \rho - 1$$
%
% $$E\left(x,t\right)= \frac{\partial\phi}{\partial x}$$
%
% where $\rho\left(x,t\right)=\int_v f(x,v,t) dv$.

format short e
addpath(genpath(pwd))

%% Step 1. Set input parameters
% pde :  A structure containing the initial condition and domain
% information. See PDE/vlasov4.m for and example. Note that this does not
% actually describe the PDE (that's done by the coefficient matricies), so
% we should probably rename this.
%
% Lev: Maximum level of the mesh in all dimensions.
%
% Deg: Degree of basis functions.
%
% TEND : End time of the simulation, i.e., run from t=0 to t=TEND in steps
% of dt.
%
% quiet : Print debugging statements or not.
%
% Compression : Choice of approach to constructing the A system matrix.

if ~exist('pde','var') || isempty(pde)
    % Equation setup
    pde = Vlasov4;
end
if ~exist('TEND','var') || isempty(TEND)
    % End time
    TEND = pde.params.TEND;
end
% if ~exist('Lev','var') || isempty(Lev)
%     % Number of levels
%     Lev = 3;
% end
if ~exist('Deg','var') || isempty(Deg)
    % Polynomial degree
    Deg = 2; % Deg = 2 Means Linear Element
end
if ~exist('quiet','var') || isempty(quiet)
    % Enable / disable print statements
    quiet = 0;
end
if ~exist('compression','var') || isempty(compression)
    % Use or not the compression reference version
    compression = 4;
end
if ~exist('implicit','var') || isempty(implicit)
    % Use or not the compression reference version
    pde.implicit = 0;
end

% Get x and v domain ranges.
Lmin = pde.params.Lmin;
Lmax = pde.params.Lmax;
Vmin = pde.params.Vmin;
Vmax = pde.params.Vmax;


% Level information.
% LevX = Lev;
% LevV = Lev;
% pde.params.LevX = LevX;
% pde.params.LevV = LevV;

% Dimensionality.
Dim = pde.Dim;
domain = pde.domain;
Lev = pde.Lev;

% DimX = 1;
% DimV = 1;
% pde.params.DimX = DimX;
% pde.params.DimV = DimV;

% Degree
pde.params.Deg = Deg;

% Time step.
dt = (domain(2,1)-domain(1,1))/2^Lev(1)/domain(2,2)/(2*Deg+1);

%% Step 1.1. Setup the multi-wavelet transform in 1D (for each dimension).
%
% Input:  Deg and Lev
%
% Output: Convert Matrix FMWT_COMP
for d = 1:Dim
    FMWT_COMP(:,:,d) = OperatorTwoScale(Deg,2^Lev(d));
    % FMWT_COMP_v = OperatorTwoScale(Deg,2^LevV);
end


%% Step 1.2. Apply the mulit-wavelet transform to the initial conditions in each dimension.
% Generate the 1D initial conditions. Input: LevX,LevV,Deg,Lmax,Vmax,pde
% Output: fval (fv and fx)--intial condition f(x,v,t=0)
%         rho--intial condition rho(x,t=0)

if ~quiet; disp('[1.2] Setting up 1D initial conditions'); end
%[fv,fx] = Intial_Con(LevX,LevV,Deg,Lmax,Vmax,pde,FMWT_COMP_x,FMWT_COMP_v);
for d = 1:Dim
    f(:,d) = forwardMWT(Lev(d),Deg,domain(1,d),domain(2,d),pde.f_0{d},pde.params);
    % fx = forwardMWT(LevX,Deg,Lmin,Lmax,pde.Fx_0,pde.params);
    % fv = forwardMWT(LevV,Deg,Vmin,Vmax,pde.Fv_0,pde.params);
end


%% Step 2. Generate Sparse-Grid (as the Hash + Connectivity tables).

%%% Construct forward and inverse hash tables.
if ~quiet; disp('[2.1] Constructing hash and inverse hash tables'); end
% [HASH,HASHInv] = HashTable(Lev,Dim);
[HASH,HASHInv] = HashTable2(Lev(1),Dim);
nHash = numel(HASHInv);

%%% Construct the connectivity.
if ~quiet; disp('[2.2] Constructing connectivity table'); end
if Dim == 1
    Con = Connect1D(Lev(1));
elseif Dim == 2
    if Lev(1)~=Lev(2)
        exit;% we need to fix the 2D connectivity
    else
        Con = Connect2D(Lev(1),HASH,HASHInv);
    end
    
end

%%% Get the multi-wavelet coefficient representation on the sparse-grid,
%%% i.e., above we transformed each of the initial condition 1D
%%% dependencies into the multi-wavelet basis. Here we combine those N 1D
%%% basis representations into the multi-D (here N=2, so 2D) initial
%%% condition, multi-wavelet representation of the initial condition
%%% specified in PDE.
if ~quiet; disp('[2.3] Calculate 2D initial condition on the sparse-grid'); end
if Dim == 1
    fval = f(:,1);
elseif Dim == 2
    fval = initial_condition_vector(f(:,1),f(:,2),Deg,Dim,HASHInv,pde);
end

clear fv fx

%% Step 3. Generate time-independent coefficient matrices
% Vlasolv Solver:
%   Operators:
%               vMassV: int_v v*l_i(v)*l_j(v)dv GradV: int_v
%               (l_i(v))'*l_j(v)dv GradX: int_x (m_i(x))'*m_j(x)dx
% Poisson Solver:
%               Operators: DelaX: int_x (m_i(x))''*m_j(x)dx
% Input:
%               LevX, LevV, k, dim, Lmax, Vmax
% Output:
%               2D Matrices--vMassV,GradV,GradX,DeltaX

%%% Build the time independent coefficient matricies.
if ~quiet; disp('[3.1] Calculate time independent matrix coefficients'); end
for d = 1:Dim
    [gMass(:,:,d),Grad(:,:,d),Delta(:,:,d)]=matrix_coeff_TI2(Lev(d),Deg,domain(1,d),domain(2,d),FMWT_COMP(:,:,d));
end
% [vMassV,GradV,GradX,DeltaX] = matrix_coeff_TI(LevX,LevV,Deg,Lmax,Vmax,...
%     FMWT_COMP_x,FMWT_COMP_v);

%%% Generate A_encode / A_data time independent data structures.
if ~quiet; disp('[3.2] Generate A_encode data structure for time independent coefficients'); end
if Dim >= 2
    if compression == 3
        A_encode=GlobalMatrixSG(gMass(:,:,2),Grad(:,:,1),HASHInv,Con,Deg);
        %     A_encode=GlobalMatrixSG(vMassV,GradX,HASHInv,Con2D,Deg);
    else
        % A_data is constructed only once per grid refinement, so can be done
        % on the host side.
        A_data = GlobalMatrixSG_SlowVersion(HASHInv,Con,Deg,compression,Dim);
    end
end

%% Step 4. Generate time-independent global matrix
% Compute the global matrix for spatial variables "x" by
%
% Poisson Solver: A_Poisson (Hash, Dim_x,k,LevX,DeltaX)
%
% Input: Hash, Dim_x,k,LevX,DeltaX,or nu, eps, CurlCurlX Output: A_Poisson
% Another Idea is to solve Poisson Equation on the finest full grid

if ~quiet; disp('[4] Construct matrix for Poisson solve'); end
% if DimX>1
%     % Construct DeltaX for DimX
% else
A_Poisson = Delta(:,:,1);
% end


%% Step 5. Time Loop
%	Step 5.1 Vlasov Equation
%       Generate time dependent coefficient matrix Generate global matrix
%       A_Vlasov(Hash,coef_mat,Dim) Apply A_Vlasov->f by RK
%	Step 5.2 Poisson Equation
%       Solve Poisson Equation sol_Poisson by A_Poisson(f) Compute
%       E=(sol_Poisson)'


% At time = 0 plotting.
if ~quiet; disp('[5.0] Plotting intial condition'); end


% Construct data for reverse MWT in 2D
% [Meval_v,v_node,Meval_x,x_node]=matrix_plot(LevX,LevV,Deg,Lmax,Vmax,...
%     FMWT_COMP_x,FMWT_COMP_v);
for d = 1:Dim
    [Meval(:,:,d),node(:,d)]=matrix_plot_1D(Lev(d),Deg,domain(1,d),domain(2,d),...
        FMWT_COMP(:,:,d));
end

if Dim == 2
    [xx,vv]=meshgrid(node(:,1),node(:,2));
    
    if ~quiet
        % Transform from wavelet space to real space
        %     tmp=Multi_2D(Meval_v,Meval_x,fval,HASHInv,Lev,Deg);
        tmp=Multi_2D(Meval(:,:,1),Meval(:,:,2),fval,HASHInv,Lev(1),Deg);
        figure(1000)
        mesh(xx,vv,reshape(tmp,Deg*2^Lev(1),Deg*2^Lev(1))','FaceColor','interp','EdgeColor','interp');
        %     axis([0 Lmax -Vmax Vmax])
        axis([domain(1,1) domain(2,1) domain(1,2) domain(2,2)])
        view(0,90)
        colorbar
    end
end% Plot initial condition


% Write the initial condition to file.
write_fval = 0;
if write_fval; write_fval_to_file(fval,Lev(1),Deg,0); end

count=1;
plotFreq = 10;

if ~quiet; disp('[7] Advancing time ...'); end
for L = 1:floor(TEND/dt)
    
    time(count) = L*dt;
    timeStr = sprintf('Step %i of %i',L,floor(TEND/dt));
    
    if ~quiet; disp(timeStr); end
    
    if pde.solvePoisson
        %%% Solve Poisson to get E (from 1-rho=1-int f dv)
        if ~quiet; disp('    [a] Solve poisson to get E'); end
        [E,u] = PoissonSolve(Lev(1),Deg,domain(1,1),domain(2,1),fval,A_Poisson,FMWT_COMP(:,:,1),domain(1,2),domain(2,2));
    end
    
    if pde.applySpecifiedE
        %%% Apply specified E
        if ~quiet; disp('    [a] Apply specified E'); end
        E = forwardMWT(Lev(1),Deg,domain(1,1),domain(2,1),pde.Ex,pde.params);
        E = E * pde.Et(time(count));
    end
    
    %%% Generate EMassX time dependent coefficient matrix.
    if ~quiet; disp('    [b] Calculate time dependent matrix coeffs'); end
    EMassX = matrix_coeff_TD(Lev(1),Deg,domain(1,1),domain(2,1),E,FMWT_COMP(:,:,1));
    
    %%% Update A_encode for time-dependent coefficient matricies.
    if ~quiet; disp('    [c] Generate A_encode for time-dependent coeffs'); end
    if compression == 3
        % broken and need to change GradV
        B_encode = GlobalMatrixSG(GradV,EMassX,HASHInv,Con2D,Deg);
        C_encode=[A_encode B_encode];
    else
        
    end
    
    %%% Advance Vlasov in time with RK3 time stepping method.
    if ~quiet; disp('    [d] RK3 time step'); end
    if compression == 3
        fval = TimeAdvance(C_encode,fval,time(count),dt,compression,Deg,pde,HASHInv);
    else
%         A_data.vMassV    = vMassV;
%         A_data.GradX     = GradX;
%         A_data.GradV     = GradV;
%         A_data.EMassX    = EMassX;
        A_data.vMassV    = gMass(:,:,2);
        A_data.GradX     = Grad(:,:,1);
        A_data.GradV     = Grad(:,:,2);
        A_data.EMassX    = EMassX;
        
        % Write the A_data structure components for use in HPC version.
        write_A_data = 0;
        if write_A_data && L==1; write_A_data_to_file(A_data,Lev(1),Deg); end
        
        fval = TimeAdvance(A_data,fval,time(count),dt,compression,Deg,pde,HASHInv);
        
    end
    
    %%% Write the present fval to file.
    if write_fval; write_fval_to_file(fval,Lev,Deg,L); end
    
    %%% Write data for FK6D test
    
    fname = ['tests/vlasov4_time_5_3/fval_',num2str(L,'%3.3i'),'.dat'];
    fd = fopen(fname,'w'); % where file.dat is the name you want to save to
    fwrite(fd,full(fval),'double'); % where U is the vector/matrix you want to store, double is the typename
    fclose(fd);
    
    %%% Plot results
    if mod(L,plotFreq)==0 && ~quiet
        if Dim == 2
            % need fixing
        figure(1000)
        
%         tmp=Multi_2D(Meval_v,Meval_x,fval,HASHInv,Lev,Deg);
        tmp=Multi_2D(Meval(:,:,1),Meval(:,:,2),fval,HASHInv,Lev(1),Deg);
        f2d = reshape(tmp,Deg*2^Lev(1),Deg*2^Lev(2))';
        mesh(xx,vv,f2d,'FaceColor','interp','EdgeColor','interp');
%         axis([0 Lmax -Vmax Vmax])
        axis([domain(1,1) domain(2,1) domain(1,2) domain(2,2)])
        view(0,90)
        colorbar
        
        title(['Time at ', timeStr])
        pause (0.01)
        end
    end
    
    %%% Check against known solution
    if pde.checkAnalytic
        
        % Check the wavelet space solution with the analytic solution
        fval_analytic = exact_solution_vector(HASHInv,pde,time(count));
        err_wavelet = sqrt(mean((fval(:) - fval_analytic(:)).^2));
        %disp(['    wavelet space absolute err : ', num2str(err_wavelet)]);
        %disp(['    wavelet space relative err : ', num2str(err_wavelet/max(abs(fval_analytic(:)))*100), ' %']);
        
        % Check the real space solution with the analytic solution
        fval_realspace = Multi_2D(Meval_v,Meval_x,fval,HASHInv,Lev,Deg);
        f2d = reshape(fval_realspace,Deg*2^LevX,Deg*2^LevV)';
        f2d_analytic = pde.ExactF(xx,vv,time(count));
        err_real = sqrt(mean((f2d(:) - f2d_analytic(:)).^2));
        %disp(['    real space absolute err : ', num2str(err_real)]);
        %disp(['    real space relative err : ', num2str(err_real/max(abs(f2d_analytic(:)))*100), ' %']);
        
        err = err_wavelet;
    end
    
    count=count+1;
    
end

end
