


%----------------------------------------------------------
% This code implements the generation of global matrix for
% 2D Poisson equation: -Delta u=f for [0,1]^2
% Maxmium Level for mesh: Np=5;
% Polynomial Degree: k=3; (quadratic element)
% Exact solution: u=sin(pi*x)*sin(pi*y)
% RHS: f=2*pi^2*sin(pi*x)*sin(pi*y)
%----------------------------------------------------------
% The following is the most inefficient implementation
% load 2D Mesh: H Index-->(nx,ix,px,ny,iy,py)
% Property for Mesh: nx+ny<=Np
%                    ix=[0:2^max(nx-1,0)-1],iy=[0:2^max(ny-1,0)-1]
%                    px=[0:k],py=[0:k]
% load 1D Mesh: M (n,i,p)--> Index for x or y component
% load Stiffness matrix S with size (DOF1DxDOF1D)
% load Mass matrix I (identity) with size (DOF1DxDOF1D)
% load RHS vector b_1D with size (DOF1Dx1)
% dof denotes the size for the 2D solver

Np=5;k=3;Dim=2;

load('Data.mat','H','S','b_1D');
I=speye(2^Np*k);

dof=size(H,1);

fprintf('dof: %i\n',dof);

A_s=sparse(dof,dof);
A_s3 = sparse(dof,dof);
A_s4 = sparse(dof,dof);
b_s=sparse(dof,1);
sol_s=sparse(dof,1);

H2 = zeros(dof,9);
H2(:,1:7)=H;

% Add 1D indices (Ix,Iy) for each element to the global H table.

for IndexI=1:dof
    
    % for x-component
    nx=H(IndexI,2);
    ix=H(IndexI,3);
    px=H(IndexI,4);
    % for y-component
    ny=H(IndexI,5);
    iy=H(IndexI,6);
    py=H(IndexI,7);
    
    Ix=map_1D(nx,ix,px,k);
    Iy=map_1D(ny,iy,py,k);
    
    H2(IndexI,8)=Ix;
    H2(IndexI,9)=Iy;
    
end

% Create global Ix and Iy variables (DOF long)

IxG = H2(:,8);
IyG = H2(:,9);

%% Assembling Matrix (for every element loop over all other elements)


% Approach 1 - Full loop

tic
fprintf('Approach 1\n');
for IndexI=1:dof
    
    % for x-component
    nx=H(IndexI,2);
    ix=H(IndexI,3);
    px=H(IndexI,4);
    % for y-component
    ny=H(IndexI,5);
    iy=H(IndexI,6);
    py=H(IndexI,7);
    
    Ix=map_1D(nx,ix,px,k);
    Iy=map_1D(ny,iy,py,k);
    
    for IndexJ=1:dof
        
        % for x-component
        mx=H(IndexJ,2);
        jx=H(IndexJ,3);
        kx=H(IndexJ,4);
        % for y-component
        my=H(IndexJ,5);
        jy=H(IndexJ,6);
        ky=H(IndexJ,7);
        
        Jx=map_1D(mx,jx,kx,k);
        Jy=map_1D(my,jy,ky,k);
        
        % pull from S and I to generate global matrix
        %  S*I+I*S
        
        val=S(Ix,Jx)*I(Iy,Jy) + I(Ix,Jx)*S(Iy,Jy);
        
        % Why is this constructed using "sparse"?
        %A_s=A_s+sparse(IndexI,IndexJ,val,dof,dof); 
        
        A_s(IndexI,IndexJ) = A_s(IndexI,IndexJ) + val;
        
    end
    
    % load 1D rhs b_1D
    val=kron(b_1D(Ix),b_1D(Iy));
    b_s(IndexI)=b_s(IndexI)+val;
end
toc

% Approach 2 - Contract inner loop to indexing with global indicies.

tic
fprintf('Approach 2\n');
for IndexI=1:dof
    
    
    A_s3(IndexI,:) = A_s3(IndexI,:) ...
        + S(IxG(IndexI),IxG).*I(IyG(IndexI),IyG) ...
        + I(IxG(IndexI),IxG).*S(IyG(IndexI),IyG);
    
end
toc

% Approach 3 - Remove loops and use only global indicies.

tic
fprintf('Approach 3\n');

A_s4 = S(IxG,IxG).*I(IyG,IyG) + I(IxG,IxG).*S(IyG,IyG);

toc


fprintf('DONE\n');

spy(A_s)

%% Solve Linear System A_s*sol_s=b_s
tic
sol_s = A_s\b_s*pi^2*2;
toc

% Load correct saved result from file and compare

load('soln.mat');
fprintf('norm(soln - correct_soln)= %f\n', norm(sol_s-sol_s2) )
fprintf('norm(A - correct_A)= %f\n', normest(A_s-A_s2) )
fprintf('norm(A3 - correct_A)= %f\n', normest(A_s3-A_s2) )
fprintf('norm(A4 - correct_A)= %f\n', normest(A_s4-A_s2) )


fprintf('END\n');





