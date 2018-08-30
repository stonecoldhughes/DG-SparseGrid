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
% clc

% Test 1
% exactf = @(x,t)(exp(t)*sin(pi*x));
% exactq = @(x,t)(-exp(t)*cos(pi*x).*(1-x.^2)*pi);
% % source = @(x,t)(exp(t)*(sin(pi*x)+2*pi*x.*cos(pi*x)+(1-x.^2)*pi^2.*sin(pi*x)));
% source = @(x)((sin(pi*x)+2*pi*x.*cos(pi*x)+(1-x.^2)*pi^2.*sin(pi*x)));
% funcCoef = @(x)( (1-x.^2) );
% funcCoef2 = @(x)( (-2*x) ); % diff(funcCoef,x)



% Test 2
exactf = @(x,t)(exp(t)*cos(pi*x));
exactq = @(x,t)(exp(t)*sin(pi*x).*(1-x.^2)*pi);
% source = @(x,t)(exp(t)*(sin(pi*x)+2*pi*x.*cos(pi*x)+(1-x.^2)*pi^2.*sin(pi*x)));
source = @(x)((cos(pi*x)-2*pi*x.*sin(pi*x)+(1-x.^2)*pi^2.*cos(pi*x)));
funcCoef = @(x)( (1-x.^2) );
funcCoef2 = @(x)( (-2*x) ); % diff(funcCoef,x)


% source = @(x)(sin(pi*x)+sin(pi*x)*pi^2);

% exactf = @(x,t)(sin(pi*x));
% exactq = @(x,t)(-10*pi*cos(pi*x));
% source = @(x)(10*sin(pi*x)*pi^2);
% funcCoef = @(x)(10);
% funcCoef2 = @(x)(x-x);

% Test 3
A = 10;
exactf = @(x,t)( (1/2)*exp(A*x)*A/sinh(A) );
exactq = @(x,t)( -(1/2)*(-x.^2+1).*exp(A*x)*A^2/sinh(A) );
% source = @(x,t)(exp(t)*(sin(pi*x)+2*pi*x.*cos(pi*x)+(1-x.^2)*pi^2.*sin(pi*x)));
source = @(x)(x-x);
funcCoef = @(x)( (1-x.^2) );
funcCoef2 = @(x)( (-2*x) ); % diff(funcCoef,x)


format short e
addpath(genpath(pwd))


Lev = 4;
Deg = 3;
num_plot = 5;




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

b = sparse(dof_1D,1);bb = sparse(dof_1D,1);
fexact = sparse(dof_1D,1);
qexact = sparse(dof_1D,1);

CFL = 0.001;
dt = CFL*h^((Deg-1)/3)/2;
% dt = CFL*h^((Deg)/3);
maxT = ceil(0.05/dt)

% Assume 
% [ I  A12]
% [A21  0 ]
% as the matrix
% e.g. we assume neuman boundary:: q=0 on boundary
% 

A21 = sparse(dof_1D,dof_1D);
% generate 1D matrix for DG
for L=0:n-1

    %---------------------------------------------
    % (funcCoef*q,d/dx p)
    %---------------------------------------------
    x0 = Lstart+L*h;
    x1 = x0+h;
    xi = quad_x*(x1-x0)/2+(x1+x0)/2;
    xmid= (x0+x1)/2;
    
    val=1/h*[Dp_val'*(quad_w.*funcCoef(xi).*p_val)]+...
             p_val'*(quad_w.*funcCoef2(xi).*p_val)/2;
         
    c = Deg*L+1:Deg*(L+1);
    
    A12 = A12 + sparse(c'*ones(1,Deg),ones(Deg,1)*c,val,dof_1D,dof_1D);
    
    val=1/h*[Dp_val'*(quad_w.*p_val)];
    A21 = A21 + sparse(c'*ones(1,Deg),ones(Deg,1)*c,val,dof_1D,dof_1D);
    
    val = sqrt(h)/2*[p_val'*(quad_w.*source(xi))]; 
    bb(c)=bb(c)+val;
    
    val = sqrt(h)/2*[p_val'*(quad_w.*exactf(xi,0))]; 
    fexact(c)=fexact(c)+val;
    
    val = sqrt(h)/2*[p_val'*(quad_w.*exactq(xi,0))]; 
    qexact(c)=qexact(c)+val;
    
    %----------------------------------------------
    % -<funcCoef*{q},p>
    %----------------------------------------------
    val=[ p_1'*funcCoef(x0)*p_2   p_1'*funcCoef(x0)*p_1,...
         -p_2'*funcCoef(x1)*p_2  -p_2'*funcCoef(x1)*p_1]/2/h;
    A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c,...
                   val(:,Deg+1:2*Deg)+val(:,2*Deg+1:3*Deg),...
                   dof_1D,dof_1D);
               
    val2 = [ p_1'*p_2   p_1'*p_1,...
            -p_2'*p_2  -p_2'*p_1]/2/h;  
    A21=A21+sparse(c'*ones(1,Deg),ones(Deg,1)*c,...
                   val2(:,Deg+1:2*Deg)+val2(:,2*Deg+1:3*Deg),...
                   dof_1D,dof_1D); 

    if L>0
        A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c-Deg,val(:,1:Deg),dof_1D,dof_1D);
        A21=A21+sparse(c'*ones(1,Deg),ones(Deg,1)*c-Deg,val2(:,1:Deg),dof_1D,dof_1D);
    elseif L == 0
%         A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c,p_1'*funcCoef(x0)*p_1/2/h,dof_1D,dof_1D);
        % periodic bc
%         A12 = A12+sparse(c'*ones(1,Deg),ones(Deg,1)*(Deg*(n-1)+1:Deg*(n)),...
%             val(:,1:Deg),dof_1D,dof_1D);
        % new implementation
        A21 = A21 +sparse(c'*ones(1,Deg),ones(Deg,1)*c,-p_1'*p_1/h/2,dof_1D,dof_1D);
        b(c) =  b(c)+p_1'*exactq(Lstart,0)/sqrt(h);

    end
    if L<n-1
        A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c+Deg,val(:,3*Deg+1:4*Deg),dof_1D,dof_1D);
        A21=A21+sparse(c'*ones(1,Deg),ones(Deg,1)*c+Deg,val2(:,3*Deg+1:4*Deg),dof_1D,dof_1D);
    elseif L == n-1

%         A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*c,-p_2'*funcCoef(x1)*p_2/2/h,dof_1D,dof_1D);
        % periodic bc
%         A12=A12+sparse(c'*ones(1,Deg),ones(Deg,1)*[1:Deg],val(:,3*Deg+1:4*Deg),dof_1D,dof_1D);
        % new implementation
        A21=A21+sparse(c'*ones(1,Deg),ones(Deg,1)*c,p_2'*p_2/h/2,dof_1D,dof_1D);
        
        b(c) =  b(c)-p_2'*exactq(Lend,0)/sqrt(h);

    end
    
    
    
    val = sqrt(h)/2*[p_val'*(quad_w.*exactf(xi,0))]; 
    f0(c) = val;
    
    
end



[quad_x,quad_w]=lgwt(num_plot,-1,1);

p_val = legendre(quad_x,Deg);
for L=0:n-1
    %---------------------------------------------
    % Generate the coefficients for DG bases
    %---------------------------------------------

    Iu = [Deg*L+1:Deg*(L+1)];

    Iv = [num_plot*L+1:num_plot*(L+1)];
%     xi=h*(quad_x/2+1/2+L);

    x0 = Lstart+L*h;
    x1 = x0+h;
    xi = quad_x*(x1-x0)/2+(x1+x0)/2;[L*h,L*h+h];%
    

    Meval(Iv,Iu)=sqrt(1/h)*p_val;
    x_node(Iv,1)=xi;

end

% checked of projection
plot(x_node,Meval*f0,'r-o',x_node,exactf(x_node,0),'b--')
hold on
plot(x_node,Meval*(A12*f0),'r-o',x_node,exactq(x_node,0),'b--')

% return
Mat = A21*A12;

% max(abs(Meval*f0))
% val = Meval*A12*f0-exactq(x_node,0);
% [norm(val) max(abs(val))]

b = b+bb;



 total_particle = 0;
 ffval = Meval*f0;

 for i = 1:num_plot
      total_particle =  total_particle+quad_w(i)*h/2*sum(ffval(i:num_plot:end));
 end
 total_particle
 
tp(1) = total_particle;

figure
for t = 1:maxT
%     t
    time = t*dt;
%     tmp = A12'*A12*f0;%dt*A12'*A12*f0+dt*b*exp(time);
%     tmp2 = A12*A12*f0;

% %     fval = f0+dt*Mat*f0+dt*b;%*exp(time);

    f1 = f0 + dt*( Mat*f0+b*exp(time-dt) );
    f2 = 3/4*f0+1/4*f1+1/4*dt*(Mat*f1+b*exp(time));
    fval = 1/3*f0+2/3*f2+2/3*dt*(Mat*f2+b*exp(time-dt/2));

%     f1 = f0 + dt*( Mat*f0+b );
%     f2 = 3/4*f0+1/4*f1+1/4*dt*(Mat*f1+b );
%     fval = 1/3*f0+2/3*f2+2/3*dt*(Mat*f2+b );
    
    f0 = fval;
    
    plot(x_node,Meval*f0,'r-o',x_node,Meval*(A12*f0),'b-<')
    pause (0.1)
    
    total_particle = 0;
    ffval = Meval*f0;
    for i = 1:num_plot
     
      total_particle =  total_particle+quad_w(i)*h/2*sum(ffval(i:num_plot:end));
    end
    tp(t+1) = total_particle;
end
% figure;plot(x,f_loc'*f0,'r-o');hold on;
% plot(x,exactf(x,time),'b--')
hold on
plot(x_node,exactf(x_node,time),'r-o')
%  max(abs(Meval*f0))
 val = Meval*f0-exactf(x_node,time);
%  [norm(val) max(abs(val))]
 
 fL2 = 0; fLinf = max(abs(val));
%  total_particle = 0;
 ffval = Meval*f0;
 for i = 1:num_plot
     fL2 = fL2 + quad_w(i)*h/2*sum(val(i:num_plot:end).^2);
%       total_particle =  total_particle+quad_w(i)*h/2*sum(ffval(i:num_plot:end));
 end
 [sqrt(fL2) fLinf]
%  total_particle
%  figure;
%  plot(x_node,Meval*A12*f0,'r-o',x_node,exactq(x_node,time),'b--')
 val = Meval*A12*f0-exactq(x_node,time);
%  [norm(val) max(abs(val))]
 
  qL2 = 0; qLinf = max(abs(val));
 for i = 1:num_plot
     qL2 = qL2 + quad_w(i)*h/2*sum(val(i:num_plot:end).^2);
     
 end
 [sqrt(qL2) qLinf]
 
%  err = f0-fexact;
%  full([norm(err) max(abs(err))])
%  
%  err = A12*f0-qexact;
%  full([norm(err) max(abs(err))])

figure;
plot(x_node,Meval*f0,'r-o',x_node,exactf(x_node,time),'r--');
hold on;
plot(x_node,Meval*A12*f0,'b-o',x_node,exactq(x_node,time),'b--');

figure;
plot(tp,'r-o')
