function [Eh,Bh,Mat,A_encode,count] = TryA_encode(Lev,Deg,Hash,InvHash,Con1D,GradX,...
                                 eps,mu,omega,dt,MaxT,...
                                 Rhs,E0,B0)
%=====================================================
% Maxwell Solver on [0,1]^3
% We use operator (du/dx,v) to construct (curl u,v)
% Note: The 3D matrix is not symmetric
%    [ 0,   A12, A13]
% As=[-A12,  0,  A23]
%    [-A13,-A23,  0 ]
%        [  0   As/(eps*nu)]
% MaxMat=[-As      0       ]
%
% MaxB=[-bs/eps,0]
% Then time stepping for
% d [Eh]        [Eh] [MaxB]
%---    =MaxMat*    -
% dt[Bh]        [Bh] [0   ]
% Input:
%       Lev,Deg,Hash,InvHash,Con1D,eps,mu,omega
%       dt,MaxT: time stepping; steps
%       GradX: Grad Matrix
%       Rhs: bs
%       E0,B0: initials for E and B
% Output:
%       Eh
%       Bh
%=====================================================
% Assembling MaxMat 

Dof_Hash = size(InvHash,2);
Dofs1 = Dof_Hash*Deg^3;
Dofs3 = 3*Dofs1;
Dofs = 2*Dofs3;

MaxMat = sparse(Dofs,Dofs);
MaxRhs = sparse (Dofs,1);
MaxSol = sparse (Dofs,1);
Mat=sparse(Dofs1,Dofs1);
count=1;
%% Assembling the global Maxwell Matrix
for i = 1:Dof_Hash
    ll = InvHash{i};

    % Lev, Cell, Index from Hash
    iLev1 = ll(1);iLev2 = ll(2);iLev3 = ll(3);
    iCel1 = ll(4);iCel2 = ll(5);iCel3 = ll(6);
    iInd1 = ll(7);iInd2 = ll(8);iInd3 = ll(9);
  
    %****************************
    % consider A12 and A21=A12'
    % I*I*T
    %****************************
    jIndex1D = Con1D{iInd3};
    
    for j = 1:size(jIndex1D,2)
        
        jLev1=iLev1;jLev2=iLev2;jLev3=ceil(log2(jIndex1D(j)));
        jCel1=iCel1;jCel2=iCel2;jCel3=max(jIndex1D(j)-1-2^(jLev3-1),0);
        
        
        if jLev1+jLev2+jLev3<=Lev
            jkey = [jLev1,jLev2,jLev3,jCel1,jCel2,jCel3];

            tmp = Hash.(sprintf('i%g_',jkey));
            tmp_mat3 = GradX(Deg*(iInd3-1)+1:Deg*iInd3,Deg*(jIndex1D(j)-1)+1:Deg*jIndex1D(j));
            tmp_mat12 = speye(Deg^2);
            tmp_mat = kron(tmp_mat12,tmp_mat3);
            
            II = Deg^3*(i-1)+1:Deg^3*i;
            JJ = Deg^3*(tmp-1)+1:Deg^3*tmp;
            A_encode{count}.IndexI=II;
            A_encode{count}.IndexJ=JJ+Dofs3+Dofs1;
            A_encode{count}.A=eye(Deg);
            A_encode{count}.B=eye(Deg);
            A_encode{count}.C=-tmp_mat3/(eps*mu);
            count=count+1;
            A_encode{count}.IndexI=II+Dofs1;
            A_encode{count}.IndexJ=JJ+Dofs3;
            A_encode{count}.A=eye(Deg);
            A_encode{count}.B=eye(Deg);
            A_encode{count}.C=tmp_mat3/(eps*mu);
            count=count+1;
            A_encode{count}.IndexI=II+Dofs3;
            A_encode{count}.IndexJ=JJ+Dofs1;
            A_encode{count}.A=eye(Deg);
            A_encode{count}.B=eye(Deg);
            A_encode{count}.C=tmp_mat3;
            count=count+1;
            A_encode{count}.IndexI=II+Dofs3+Dofs1;
            A_encode{count}.IndexJ=JJ;
            A_encode{count}.A=eye(Deg);
            A_encode{count}.B=eye(Deg);
            A_encode{count}.C=-tmp_mat3;
            count=count+1;
%             [JJ,II] = meshgrid(JJ,II);
%             Iv = [II,II+Dofs1,II+Dofs3,II+Dofs3+Dofs1];
%             Iu = [JJ+Dofs3+Dofs1,JJ+Dofs3,JJ+Dofs1,JJ];
%             Aij = -[tmp_mat/(eps*mu),...
%                   -tmp_mat/(eps*mu),...
%                   -tmp_mat,...
%                    tmp_mat];
% 
%             MaxMat = MaxMat+sparse(Iv,Iu,Aij,Dofs,Dofs);
%             Mat=Mat-sparse(II,JJ,tmp_mat,Dofs1,Dofs1);
        end
        
        
    end

    %****************************
    % consider A13 and A31=A13'
    % I*T*I
    %****************************
    jIndex1D = Con1D{iInd2};

    for j = 1:size(jIndex1D,2)
        
        jLev1=iLev1;jLev2=ceil(log2(jIndex1D(j)));jLev3=iLev3;
        jCel1=iCel1;jCel2=max(jIndex1D(j)-1-2^(jLev2-1),0);jCel3=iCel3;
        
        if jLev1+jLev2+jLev3<=Lev
            jkey = [jLev1,jLev2,jLev3,jCel1,jCel2,jCel3];
            % tmp is not contigent
            tmp = Hash.(sprintf('i%g_',jkey));
            tmp_mat2 = GradX(Deg*(iInd2-1)+1:Deg*iInd2,Deg*(jIndex1D(j)-1)+1:Deg*jIndex1D(j));
            tmp_mat1 = eye(Deg);tmp_mat3=eye(Deg);
            tmp_mat = kron(kron(tmp_mat1,tmp_mat2),tmp_mat3);
            
            II = Deg^3*(i-1)+1:Deg^3*i;
            JJ = Deg^3*(tmp-1)+1:Deg^3*tmp;
            
            A_encode{count}.IndexI=II;
            A_encode{count}.IndexJ=JJ+Dofs3+2*Dofs1;
            A_encode{count}.A=eye(Deg);
            A_encode{count}.B=tmp_mat2/(eps*mu);
            A_encode{count}.C=eye(Deg);
            count=count+1;
            A_encode{count}.IndexI=II+2*Dofs1;
            A_encode{count}.IndexJ=JJ+Dofs3;
            A_encode{count}.A=eye(Deg);
            A_encode{count}.B=-tmp_mat2/(eps*mu);
            A_encode{count}.C=eye(Deg);
            count=count+1;
            A_encode{count}.IndexI=II+Dofs3;
            A_encode{count}.IndexJ=JJ+2*Dofs1;
            A_encode{count}.A=eye(Deg);
            A_encode{count}.B=-tmp_mat2;
            A_encode{count}.C=eye(Deg);
            count=count+1;
            A_encode{count}.IndexI=II+Dofs3+2*Dofs1;
            A_encode{count}.IndexJ=JJ;
            A_encode{count}.A=eye(Deg);
            A_encode{count}.B=tmp_mat2;
            A_encode{count}.C=eye(Deg);
            count=count+1;
%             [JJ,II] = meshgrid(JJ,II);
% 
%             Iv = [II,II+2*Dofs1,II+Dofs3,II+Dofs3+2*Dofs1];
%             Iu = [JJ+Dofs3+2*Dofs1,JJ+Dofs3,JJ+2*Dofs1,JJ];
%             Aij =  [tmp_mat/(eps*mu),...
%                    -tmp_mat/(eps*mu),...
%                    -tmp_mat,...
%                     tmp_mat];
% 
%             MaxMat = MaxMat+sparse(Iv,Iu,Aij,Dofs,Dofs);
%             
            
        end
    end
    
    %****************************
    % consider A23 and A32=A23'
    % T*I*I
    %****************************
    jIndex1D = Con1D{iInd1};

    for j = 1:size(jIndex1D,2)
        
        jLev1=ceil(log2(jIndex1D(j)));jLev2=iLev2;jLev3=iLev3;
        jCel1=max(jIndex1D(j)-1-2^(jLev1-1),0);jCel2=iCel2;jCel3=iCel3;
        
        if jLev1+jLev2+jLev3<=Lev
            jkey = [jLev1,jLev2,jLev3,jCel1,jCel2,jCel3];
            % tmp is not contigent
            tmp = Hash.(sprintf('i%g_',jkey));
            tmp_mat1 = GradX(Deg*(iInd1-1)+1:Deg*iInd1,Deg*(jIndex1D(j)-1)+1:Deg*jIndex1D(j));
            tmp_mat23 = speye(Deg^2);
            tmp_mat = kron(tmp_mat1,tmp_mat23);
            
            II = Deg^3*(i-1)+1:Deg^3*i;
            JJ = Deg^3*(tmp-1)+1:Deg^3*tmp;
            
            A_encode{count}.IndexI=II+Dofs1;
            A_encode{count}.IndexJ=JJ+Dofs3+2*Dofs1;
            A_encode{count}.A=-tmp_mat1/(eps*mu);
            A_encode{count}.B=eye(Deg);
            A_encode{count}.C=eye(Deg);
            count=count+1;
            A_encode{count}.IndexI=II+2*Dofs1;
            A_encode{count}.IndexJ=JJ+Dofs3+Dofs1;
            A_encode{count}.A=tmp_mat1/(eps*mu);
            A_encode{count}.B=eye(Deg);
            A_encode{count}.C=eye(Deg);
            count=count+1;
            A_encode{count}.IndexI=II+Dofs3+Dofs1;
            A_encode{count}.IndexJ=JJ+2*Dofs1;
            A_encode{count}.A=tmp_mat1;
            A_encode{count}.B=eye(Deg);
            A_encode{count}.C=eye(Deg);
            count=count+1;
            A_encode{count}.IndexI=II+Dofs3+2*Dofs1;
            A_encode{count}.IndexJ=JJ+Dofs1;
            A_encode{count}.A=-tmp_mat1;
            A_encode{count}.B=eye(Deg);
            A_encode{count}.C=eye(Deg);
            count=count+1;
%             [JJ,II] = meshgrid(JJ,II);
% 
%             Iv = [II+Dofs1,II+2*Dofs1,II+Dofs3+Dofs1,II+Dofs3+2*Dofs1];
%             Iu = [JJ+Dofs3+2*Dofs1,JJ+Dofs3+Dofs1,JJ+2*Dofs1,JJ+Dofs1];
%             Aij = -[tmp_mat/(eps*mu),...
%                   -tmp_mat/(eps*mu),...
%                   -tmp_mat,...
%                    tmp_mat];
% 
%             MaxMat = MaxMat+sparse(Iv,Iu,Aij,Dofs,Dofs);
            
        end
    end
    
    
end

%% Time advance for solving Maxwell equation
MaxRhs = -[Rhs/(eps);zeros(Dofs3,1)];
MaxSol = [E0;B0];

time=0;
for T=1:MaxT
    time=time+dt;

    bbb=MaxRhs*sin(omega*time);

    MaxSol = TimeAdvance(A_encode,MaxSol,dt,bbb);
    
end

Eh = MaxSol(1:Dofs3);
Bh = MaxSol(Dofs3+1:end);