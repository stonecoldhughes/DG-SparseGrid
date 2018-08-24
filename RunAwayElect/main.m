% main code for LDG Poisson equation

% df/dt = d/dx(1-x^2)df/dx+source
% Method 1. LDG
% [A11 A12]
% [A21 A22]
% Here A11 = I, A12 = -(d/dx sqrt(1-x^2)q,p)
% A21 = -(sqrt(1-x^2)df/dx,w), A22 = (q,w) = I
% A12 = (sqrt(1-x^2)q,dp/dx)-<sqrt(1-x^2)\hat{q},p>
% A21 = (sqrt(1-x^2)f,dw/dx)-<sqrt(1-x^2)\hat{f},w>
% try with central flux
clear all
close all
clc

exactf = @(x,t)(exp(t)*sin(pi*x));
% source = @(x,t)(exp(t)*(sin(pi*x)+2*pi*x.*cos(pi*x)+(1-x.^2)*pi^2.*sin(pi*x)));
% source = @(x)((sin(pi*x)+2*pi*x.*cos(pi*x)+(1-x.^2)*pi^2.*sin(pi*x)));
% funcCoef = @(x)(sqrt(1-x.^2));

% source = @(x)(sin(pi*x)+sin(pi*x)*pi^2);

exactf = @(x,t)(sin(pi*x));
source = @(x)(sin(pi*x)*pi^2);
funcCoef = @(x)(1);


format short e
addpath(genpath(pwd))

Lev = 7;
Deg = 1;



Lstart = -1;
Lend = 1;
Lmax = Lend-Lstart;


%--Quadrature
quad_num=10;
%---------------

% compute the trace values
p_1 = legendre(-1,Deg);
p_2 = legendre(1,Deg);

[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);
Dp_val = dlegendre(quad_x,Deg);

%---------------------------
% Jacobi of variable x and v
% Define Matrices
%---------------------------
n=2^(Lev);h=Lmax/n;
Jacobi=h;
dof_1D=Deg*n;
A12 = sparse(dof_1D,dof_1D);
f0 = sparse(dof_1D,1);

b = sparse(dof_1D,1);

CFL = 0.01;
dt = CFL*h;


% generate 1D matrix for DG
for L=0:n-1

    %---------------------------------------------
    % (funcCoef*q,d/dx p)
    %---------------------------------------------
    x0 = Lstart+L*h;
    x1 = x0+h;
    xi = quad_x*(x1-x0)/2+(x1+x0)/2;
    
    val=1/h*[Dp_val'*(quad_w.*funcCoef(xi).*p_val)];
    c = Deg*L+1:Deg*(L+1);
    
    A12 = A12 + sparse(c'*ones(1,Deg),ones(Deg,1)*c,val,dof_1D,dof_1D);
    
    %----------------------------------------------
    % -<funcCoef*{q},p>
    %----------------------------------------------
    val=[ p_1'*funcCoef(x0)*p_2    p_1'*funcCoef(x0)*p_1,...
             -p_2'*funcCoef(x1)*p_2  -p_2'*funcCoef(x1)*p_1]/2;
    A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c,...
                   val(:,Deg+1:2*Deg)+val(:,2*Deg+1:3*Deg),...
                   dof_1D,dof_1D);

    if L>0
        A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c-Deg,val(:,1:Deg),dof_1D,dof_1D);
    elseif L == 0
        %A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c,p_1'*funcCoef(x0)*p_1/2,dof_1D,dof_1D);
        % periodic bc
        A12 = A12+sparse(c'*ones(1,Deg),ones(Deg,1)*(Deg*(n-1)+1:Deg*(n)),...
            val(:,1:Deg),dof_1D,dof_1D);
    end
    if L<n-1
        A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c+Deg,val(:,3*Deg+1:4*Deg),dof_1D,dof_1D);
    elseif L == n-1
        %A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c,p_2'*funcCoef(x1)*p_2/2,dof_1D,dof_1D);
        % periodic bc
        A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*[1:Deg],val(:,3*Deg+1:4*Deg),dof_1D,dof_1D);
    end
    
    val = sqrt(h)/2*[p_val'*(quad_w.*source(xi))]; 
    b(c)=val;
    
    val = sqrt(h)/2*[p_val'*(quad_w.*exactf(xi,0))]; 
    f0(c) = val;
    
    
end

% x = [Lstart:0.01:Lend];
% [f_loc] = EvalWavPoint4(Lstart,Lend,Lev,Deg,x,2);
% plot(f0)
% return
[quad_x,quad_w]=lgwt(Deg,-1,1);
p_val = legendre(quad_x,Deg);
for L=0:n-1
    %---------------------------------------------
    % Generate the coefficients for DG bases
    %---------------------------------------------
    Iu=[Deg*L+1:Deg*(L+1)];
%     xi=h*(quad_x/2+1/2+L);
    x0 = Lstart+L*h;
    x1 = x0+h;
    xi = quad_x*(x1-x0)/2+(x1+x0)/2;
    
    Meval(Iu,Iu)=sqrt(1/h)*p_val;
    x_node(Deg*L+1:Deg*L+Deg,1)=xi;

end

% checked of projection
plot(x_node,Meval*f0,'r-o',x_node,exactf(x_node,0),'b--')
% hold on
figure
plot(x_node,source(x_node),'r-o',x_node,Meval*b,'b--')
% return
Mat = A12*A12;
figure
for t = 1:100
    time = t*dt;
%     tmp = A12'*A12*f0;%dt*A12'*A12*f0+dt*b*exp(time);
%     tmp2 = A12*A12*f0;

    fval = f0+dt*Mat*f0+dt*b;%*exp(time);
    
%     f1 = f0 + dt*( Mat*f0+b*exp(time-dt) );
%     f2 = 3/4*f0+1/4*f1+1/4*dt*(Mat*f1+b*exp(time));
%     fval = 1/3*f0+2/3*f2+2/3*dt*(Mat*f2+b*exp(time-dt/2));

%     f1 = f0 + dt*( Mat*f0+b );
%     f2 = 3/4*f0+1/4*f1+1/4*dt*(Mat*f1+b );
%     fval = 1/3*f0+2/3*f2+2/3*dt*(Mat*f2+b );
    
    f0 = fval;
    
    plot(x_node,Meval*f0)
%    plot(x_node,Meval*(tmp),'r-o',x_node,Meval*(tmp2),'b-<')
    pause(0.1)
end
% figure;plot(x,f_loc'*f0,'r-o');hold on;
% plot(x,exactf(x,time),'b--')
hold on
plot(x_node,exactf(x_node,time),'r-o')


