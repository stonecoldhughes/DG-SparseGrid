function bc = ComputeBC(pde,t,Lev,Deg,xMin,xMax,Fun,LorR)
% function bc = ComputeBC(Lev,Deg,xMin,xMax,Fun,time,LorR)
% function ComputeBC to compute the bc term
% This is the evaluation for two points on 1D
% Func*v|_xMin and Func*v|_xMax
%-----------------------------------------------------

p = pde.params;

L = xMax-xMin;
Tol_Cel_Num = 2^(Lev);
h = L  / Tol_Cel_Num;
DoF = Deg * Tol_Cel_Num;

bc = sparse(DoF,1);

p_L = legendre(-1,Deg) * 1/sqrt(h);
p_R = legendre(+1,Deg) * 1/sqrt(h);

if strcmp(LorR,'L')
    
    WorkCel = 0;
    c = [1:Deg];
    IntVal =  p_L'*Fun(xMin,p,t) ;
%     IntVal =  p_L'*(Fun(xMin,time)) ;

    bc(c) = - IntVal;
    
else
    
    WorkCel = Tol_Cel_Num - 1;
    c = Deg*WorkCel+[1:Deg];
    IntVal =  p_R'*Fun(xMax,p,t);
%     IntVal =  p_R'*(Fun(xMax,time));
    bc(c) = IntVal;
    
end

end
