function Mat = MatrixGrad(Lev,Deg,LInt,LEnd,FunCoef)
%function Mat to compute the Grad operator
%-----------------------------------------------------
L = LEnd-LInt;
Tol_Cel_Num = 2^(Lev);
h = L  / Tol_Cel_Num;
DoF = Deg * Tol_Cel_Num;

Mat = sparse(DoF,DoF);

quad_num = 10;
FluxVal = 1;
% compute the trace values
p_L = legendre(-1,Deg) * 1/sqrt(h);
p_R = legendre(+1,Deg) * 1/sqrt(h);

%%
%  Get the basis functions and derivatives for all k
%  p_val(:,:) is quad_num by deg
%  Dp_val(:,:) is quad_num by deg
[quad_x,quad_w] = lgwt(quad_num,-1,1);
p_val  = legendre(quad_x,Deg)  * 1/sqrt(h);
Dp_val = dlegendre(quad_x,Deg) * 1/sqrt(h) * 2/h;

for WorkCel = 0 : Tol_Cel_Num - 1
    %---------------------------------------------
    % (funcCoef*q,d/dx p)
    %---------------------------------------------
    c = Deg*WorkCel+[1:Deg];
    
    xL = LInt + WorkCel*h;
    xR = xL + h;
    PhyQuad = quad_x*(xR-xL)/2+(xR+xL)/2;
    
    IntVal=[Dp_val'*(quad_w.*FunCoef(PhyQuad).*p_val)];
    Mat = Mat + sparse(c'*ones(1,Deg),ones(Deg,1)*c,IntVal,DoF,DoF);
    %----------------------------------------------
    % -<funcCoef*{q},p>
    %----------------------------------------------
%     TraVal = [...
%               -p_L' * FunCoef(xL) * ( p_R/2 - FluxVal/2*p_R),...
%               -p_L' * FunCoef(xL) * ( p_L/2 + FluxVal/2*p_L),... % xL
%                p_R' * FunCoef(xR) * ( p_R/2 - FluxVal/2*p_R),...
%                p_R' * FunCoef(xR) * ( p_L/2 + FluxVal/2*p_L),... % xR
%             ];
    TraVal = [...
              -p_L' *  ( FunCoef(xL) * p_R/2 - FluxVal/2*abs(FunCoef(xL))*p_R),...
              -p_L' *  ( FunCoef(xL) * p_L/2 + FluxVal/2*abs(FunCoef(xL))*p_L),... % xL
               p_R' *  ( FunCoef(xR) * p_R/2 - FluxVal/2*abs(FunCoef(xR))*p_R),...
               p_R' *  ( FunCoef(xR) * p_L/2 + FluxVal/2*abs(FunCoef(xR))*p_L),... % xR
            ];

    % Adding trace value to matrix
    RowInd = [c'*ones(1,Deg) c'*ones(1,Deg) c'*ones(1,Deg) c'*ones(1,Deg)];
    ColInd = [ones(Deg,1)*(c-Deg),ones(Deg,1)*c,ones(Deg,1)*c,ones(Deg,1)*(c+Deg)];
    
    if WorkCel == 0
        Iu = RowInd(:,Deg+1:end);
        Iv = ColInd(:,Deg+1:end);
        Val = TraVal(:,Deg+1:end);
    elseif WorkCel == Tol_Cel_Num - 1
        Iu = RowInd(:,1:3*Deg);
        Iv = ColInd(:,1:3*Deg);
        Val = TraVal(:,1:3*Deg);
    else
        Iu = RowInd;
        Iv = ColInd;
        Val = TraVal;
    end
    Mat = Mat + sparse(Iu,Iv,Val,DoF,DoF);
    
%     Mat = Mat + sparse(RowInd,ColInd,TraVal,DoF,DoF);
    
    
end

figure;spy(Mat)

end