clear
clc
close all

sigma = 0.1;
x0 = 0.3;
f0 = @(x)(exp(-(x-x0).^2/sigma^2)+exp(-(x+x0).^2/sigma^2));
exactf = @(x,t)(f0(x));

format short e
addpath(genpath(pwd))


Lev = 10;
Deg = 1;
num_plot = 3;

PlotType = 1;
AdapTest = 1;

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

tol_cel_num=2^(Lev);h=Lmax/tol_cel_num;
dof_1D=Deg*tol_cel_num;

for L=0:tol_cel_num-1
    
    %---------------------------------------------
    % (funcCoef*q,d/dx p)
    %---------------------------------------------
    x0 = Lstart+L*h;
    x1 = x0+h;
    xi = quad_x*(x1-x0)/2+(x1+x0)/2;
    xmid= (x0+x1)/2;
    c = Deg*L+1:Deg*(L+1);
    
    val = sqrt(h)/2*[p_val'*(quad_w.*exactf(xi,0))];
    fval(c,1) = val;
    
    
end

% this is hard for p adaptive
% % h = 2/4;
% % Deg_loc=[1 2 2 1];
% % count = 0;
% % for L=0:3
% %     
% %     %---------------------------------------------
% %     % (funcCoef*q,d/dx p)
% %     %---------------------------------------------
% %     x0 = Lstart+L*h;
% %     x1 = x0+h;
% %     xi = quad_x*(x1-x0)/2+(x1+x0)/2;
% %     xmid= (x0+x1)/2;
% %     c = count+1:count+Deg_loc(L+1);%Deg_loc(L+1)*L+1:Deg_loc(L+1)*(L+1);
% %     count = count+Deg_loc(L+1);
% %     p_val = legendre(quad_x,Deg_loc(L+1));
% %     Dp_val = dlegendre(quad_x,Deg_loc(L+1));
% % 
% %     val = sqrt(h)/2*[p_val'*(quad_w.*exactf(xi,0))];
% %     fval_p(c,1) = val;
% %     
% % 
% % end
% % 
% % return
[quad_x,quad_w]=lgwt(num_plot,-1,1);

p_val = legendre(quad_x,Deg);
for L=0:tol_cel_num-1
    %---------------------------------------------
    % Generate the coefficients for DG bases
    %---------------------------------------------
    
    Iu = [Deg*L+1:Deg*(L+1)];
    
    Iv = [num_plot*L+1:num_plot*(L+1)];
    
    x0 = Lstart+L*h;
    x1 = x0+h;
    
    if PlotType == 1
        xi = [x0,quad_x(2:end-1)'*(x1-x0)/2+(x1+x0)/2,x1];
    elseif PlotType == 0
        xi = quad_x*(x1-x0)/2+(x1+x0)/2;
    end
    
    
    Meval(Iv,Iu)=sqrt(1/h)*p_val;
    x_node(Iv,1)=xi;
    
end

FMWT = OperatorTwoScale(Deg,2^Lev);
%         Mat = FMWT*Mat*FMWT';
%         b = FMWT*b;
Meval = Meval*FMWT';
%
fval = FMWT*fval;
% return
Lev = 3;
clear fcell
ffval = fval(1:2^Lev);
for i = 0:Lev
    if i == 0
        startP = 1;
    else
        startP = Deg*(2^max(i-1,0))+1;
    end
    %     endP = Deg*2^i;
    %     fcell{i+1} = ffval(startP:endP);
    for cell = 0:2^max(i-1,0)-1
        
        fcell{i+1}{cell+1} = {ffval(startP+cell)};
    end
    
    
end
% return
%
%         save(['Mat_Lev',num2str(Lev),'_Deg',num2str(Deg),'.mat']);

% checked of projection
plot(x_node(1:2^4:end),Meval(1:2^4:end,1:2^Lev)*ffval,'r-o',x_node,exactf(x_node,0),'b--','LineWidth',2)
legend({'solution f','flux x(1-x^2)f'})
title(['time at ',num2str(0)])
% return
[Hash,IHash,FineIndex] = HashTable1D(Lev);
subplot(1,2,1)
plotgrid(IHash,Lstart,Lend,0,'r-o');hold on;
subplot(1,2,2)
Id = [];
for i = 1:size(IHash,2)
Id = [Id;IHash{i}(3)];
end
plot(x_node(1:2^4:end),Meval(1:2^4:end,Id)*ffval,'r-o',x_node,exactf(x_node,0),'b--','LineWidth',2)
% plot(x_node,Meval(:,Id)*ffval-exactf(x_node,0),'b--','LineWidth',2)
hold on;

DoFs = Deg*2^Lev;
MMeval = Meval(:,1:DoFs);


eps = 1e-2;%1.1636e-2;
% from the finest grid
FineStart = DoFs-2^(Lev-1)+1;
oldIHash = IHash;
for iter = 1:10
    count = size(IHash,2);
    pp = size(oldIHash,2);
    for index = FineStart:DoFs
        
        
        ll = IHash{index};
        lev_loc = ll(1);
        cel_loc = ll(2);
        
        f_loc = ffval(index);
%         [lev_loc,cel_loc]
        key = [lev_loc-1,floor(cel_loc/2)];
        mom_index = Hash.(sprintf('i%g_',key));
        f_mom = ffval(mom_index);
        if abs(f_loc)> eps
%         if abs(f_loc)>abs(f_mom)/2
%             'refine'
            key = [lev_loc+1,2*cel_loc];
            if isfield(Hash,sprintf('i%g_',key)) == 0
            count = count+1;
            Hash = setfield(Hash,sprintf('i%g_',key),count);
            IHash{count} = [key,2^(lev_loc)+2*cel_loc+1];
            ffval(count)=fval(2^(lev_loc)+2*cel_loc+1);
            end
            
            key = [lev_loc+1,2*cel_loc+1];
            if isfield(Hash,sprintf('i%g_',key)) == 0
            count = count+1;
            Hash = setfield(Hash,sprintf('i%g_',key),count);
            IHash{count} = [key,2^(lev_loc)+2*cel_loc+2];
            ffval(count)=fval(2^(lev_loc)+2*cel_loc+2);
            end
            
            
            
            
            
%         else
%             'exit refine'
        end
        
    end
    subplot(1,2,1)
    plotgrid(IHash,Lstart,Lend,iter,'b-o');
    oldIHash = IHash;
    FineStart = pp+1;
    DoFs = size(IHash,2);
    
    Id = [];
for i = 1:size(IHash,2)
Id = [Id;IHash{i}(3)];
end
subplot(1,2,2)
plot(x_node(1:2^4:end),Meval(1:2^4:end,Id)*ffval+iter,'r-o',x_node,exactf(x_node,0)+iter,'b--','LineWidth',2)
% plot(x_node,Meval(:,Id)*ffval-exactf(x_node,0)+iter,'b--','LineWidth',2)


    [iter FineStart DoFs]
    if FineStart >= DoFs
        break
    end
end
% axis([-1 1 0 10])

ii = setdiff([1:size(fval,1)],Id);
fnot = fval(ii);
max(abs(fnot))