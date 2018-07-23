function [B,E1,E2] = SemiTime(Mat,dt,B,E1,E2,b1,b2) %Backward Euler
%-------------------------------------------------
% Time Advance Method
% Input: Matrix:: A
%        Vector:: f
%        Time Step:: dt
% Output: Vector:: f
%-------------------------------------------------
[B,E1,E2] = Semi(Mat,dt,B,E1,E2,b1,b2);

end

function [B,E1,E2] = Semi(Mat,dt,B,E1,E2,b1,b2)
%----------------------------------
% Semi-Implicit Method
%----------------------------------
B=Mat*E1*dt+B;
E1=Mat*B*dt+dt*b1+E1;
E2=dt*b2+E2;

%     sol_1=sol_n+dt*(A_s*sol_n+b_s*sin(pde.w*time));
%     sol_2=3/4*sol_n+1/4*sol_1+1/4*dt*(A_s*sol_1+b_s*sin(pde.w*time));
%     sol_n=1/3*sol_n+2/3*sol_2+2/3*dt*(A_s*sol_2+b_s*sin(pde.w*time));
end

function ftmp = ApplyA(A,f)
%-----------------------------------
% Multiply Matrix A by Vector f
%-----------------------------------
dof = size(f,1);
ftmp=sparse(dof,1);

ftmp = A*f;
% for i=1:size(A,2)
% tmpA=A{i}.A1;
% tmpB=A{i}.A2;
% IndexI=A{i}.IndexI;
% IndexJ=A{i}.IndexJ;
% ftmp(IndexI)=ftmp(IndexI)+kron_mult2(tmpA,tmpB,f(IndexJ));
% end

end






