function [Eh,Bh] = MaxwellSolver3(Lev,Deg,Hash,InvHash,Con1D,Grad1X,Grad2X,...
                                 eps,mu,omega,dt,MaxT,...
                                 Rhs,E0,B0)
%=====================================================
% Maxwell Solver on [0,1]^3
% We use operator (du/dx,v) to construct (curl u,v)
% Note: The 3D matrix is not symmetric
%    [ 0,   A12, A13]
% As=[-A12,  0,  A23]
%    [-A13,-A23,  0 ]
%        [ Bs/ (eps*mu)  As/(eps*mu)]
% MaxMat=[-As                 Bs    ]
%    [ B11,  0,   0 ]
% Bs=[ 0,   B22,  0 ]
%    [ 0,    0,  B33]
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
            tmp_mat3 = Grad1X(Deg*(iInd3-1)+1:Deg*iInd3,Deg*(jIndex1D(j)-1)+1:Deg*jIndex1D(j));
            tmp_mat12 = speye(Deg^2);
            tmp_mat = kron(tmp_mat12,tmp_mat3); %"-"sign
            
            
            II = Deg^3*(i-1)+1:Deg^3*i;
            JJ = Deg^3*(tmp-1)+1:Deg^3*tmp;
            [II,JJ] = meshgrid(JJ,II);
            
            Iu1 = [II,II+Dofs1,II+Dofs3,II+Dofs3+Dofs1];
            Iv1 = [JJ+Dofs3+Dofs1,JJ+Dofs3,JJ+Dofs1,JJ];
            Aij = [tmp_mat/(eps*mu),...
                  -tmp_mat/(eps*mu),...
                  -tmp_mat,...
                   tmp_mat];
            MaxMat = MaxMat+sparse(Iu1,Iv1,Aij,Dofs,Dofs);
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
            tmp_mat2 = Grad1X(Deg*(iInd2-1)+1:Deg*iInd2,Deg*(jIndex1D(j)-1)+1:Deg*jIndex1D(j));
            tmp_mat1 = eye(Deg);tmp_mat3=eye(Deg);
            tmp_mat = kron(kron(tmp_mat1,tmp_mat2),tmp_mat3);
            
            
            
            II = Deg^3*(i-1)+1:Deg^3*i;
            JJ = Deg^3*(tmp-1)+1:Deg^3*tmp;
            [II,JJ] = meshgrid(JJ,II);

            Iu1 = [II,II+2*Dofs1,II+Dofs3,II+Dofs3+2*Dofs1];
            Iv1 = [JJ+Dofs3+2*Dofs1,JJ+Dofs3,JJ+2*Dofs1,JJ];
            Aij = -[tmp_mat/(eps*mu),...
                   -tmp_mat/(eps*mu),...
                   -tmp_mat,...
                    tmp_mat];
            MaxMat = MaxMat+sparse(Iu1,Iv1,Aij,Dofs,Dofs);
            
            
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
            tmp_mat1 = Grad1X(Deg*(iInd1-1)+1:Deg*iInd1,Deg*(jIndex1D(j)-1)+1:Deg*jIndex1D(j));
            tmp_mat23 = speye(Deg^2);
            tmp_mat = kron(tmp_mat1,tmp_mat23); %"-"sign
         
            
            II = Deg^3*(i-1)+1:Deg^3*i;
            JJ = Deg^3*(tmp-1)+1:Deg^3*tmp;
            [II,JJ] = meshgrid(JJ,II);

            Iu1 = [II+Dofs1,II+2*Dofs1,II+Dofs3+Dofs1,II+Dofs3+2*Dofs1];
            Iv1 = [JJ+Dofs3+2*Dofs1,JJ+Dofs3+Dofs1,JJ+2*Dofs1,JJ+Dofs1];
            Aij = [tmp_mat/(eps*mu),...
                  -tmp_mat/(eps*mu),...
                  -tmp_mat,...
                   tmp_mat];
            MaxMat = MaxMat+sparse(Iu1,Iv1,Aij,Dofs,Dofs);
            
        end
    end
    
   jIndex1 = Con1D{iInd1};
   jIndex2 = Con1D{iInd2};
   jIndex3 = Con1D{iInd3};
    
   for j1 = 1:size(jIndex1,2)
       for j2 = 1:size(jIndex2,2)
           for j3 = 1:size(jIndex3,2)
               jLev1=ceil(log2(jIndex1(j1)));jLev2=ceil(log2(jIndex2(j2)));jLev3=ceil(log2(jIndex3(j3)));
               jCel1=max(jIndex1(j1)-1-2^(jLev1-1),0);jCel2=max(jIndex2(j2)-1-2^(jLev2-1),0);jCel3=max(jIndex3(j3)-1-2^(jLev3-1),0);
               
               
           if jLev1+jLev2+jLev3<=Lev
                  jkey = [jLev1,jLev2,jLev3,jCel1,jCel2,jCel3];
                  tmp = Hash.(sprintf('i%g_',jkey));
                  up_mat1 = Grad2X(Deg*(iInd1-1)+1:Deg*iInd1,Deg*(jIndex1(j1)-1)+1:Deg*jIndex1(j1));
                  up_mat2 = Grad2X(Deg*(iInd2-1)+1:Deg*iInd2,Deg*(jIndex2(j2)-1)+1:Deg*jIndex2(j2));
                  up_mat3 = Grad2X(Deg*(iInd3-1)+1:Deg*iInd3,Deg*(jIndex3(j3)-1)+1:Deg*jIndex3(j3));
                  up=kron(kron(up_mat1,up_mat2),up_mat3);
                  
                  II = Deg^3*(i-1)+1:Deg^3*i;
                  JJ = Deg^3*(tmp-1)+1:Deg^3*tmp;
                  [II,JJ] = meshgrid(JJ,II);
                  
                  Iu=[II,II+Dofs1,II+2*Dofs1,II+Dofs3,II+Dofs3+Dofs1,II+Dofs3+2*Dofs1];
                  Iv=[JJ,JJ+Dofs1,JJ+2*Dofs1,JJ+Dofs3,JJ+Dofs3+Dofs1,JJ+Dofs3+2*Dofs1];
                  Bij=[up,up,up,up,up,up]; 
                  MaxMat = MaxMat+sparse(Iu,Iv,Bij,Dofs,Dofs);
           end
           end
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

    MaxSol = TimeAdvance(MaxMat,MaxSol,dt,bbb);
    
end

Eh = MaxSol(1:Dofs3);
Bh = MaxSol(Dofs3+1:end);
