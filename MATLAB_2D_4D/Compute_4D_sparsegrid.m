% Computing on the sparse grid for Poisson Eq
% 4-dimensional calculation
%------------------------------------------------
% A_s: Matrix
% b_s: RHS
% sol_s: Solution
% uu_s: Interpolation
% Loop sum_level
%------------------------------------------------

dim=4;

% Key1DMesh
% ['Generate 1D keys']
% tic
nx=[];px=[];kx=[];
for Lx=0:n
    for Px=0:2^max(0,Lx-1)-1
        for Kx=1:k
            nx=[nx;Lx];
            px=[px;Px];
            kx=[kx;Kx];
            
        end
    end
end
%
Key1dMesh=[nx,px,kx];
% toc

['Generate ',num2str(dim),'-D Hash Table']
% tic
Hash =Mesh_interface;
Hash.init (dim,n,k-1,k^dim);
% toc
dof_sparse = double(Hash.size);

% sort the index
All_index=allv(Hash);
All_index=sort(All_index);
Real_index=All_index'+1;

A_s = sparse(dof_sparse,dof_sparse);
b_s = sparse(dof_sparse,1);
sol_s = sparse(dof_sparse,1);
uu_s=sparse(dof_sparse,1);

total_kron_flops = 0;
% Method
for sum_level=0:n
    for i1_level=0:sum_level
        for i2_level=0:sum_level-i1_level
            for i3_level=0:sum_level-i1_level-i2_level
                i4_level=sum_level-i1_level-i2_level-i3_level;
                
                I1=Index_1D(k,i1_level);
                I2=Index_1D(k,i2_level);
                I3=Index_1D(k,i3_level);
                I4=Index_1D(k,i4_level);
                
                
                
                key_i=GenerateKey(I1,I2,I3,I4,Key1dMesh);
                Index_I=findvs(Hash,key_i);
                
                % Term 1: S*I*I*I+I*S*I*I+I*I*S*I+I*I*I*S at the diagonal entry
                tmp=kron(kron(kron(Stiff_1D(I1,I1),M_mass(I2,I2)),M_mass(I3,I3)),M_mass(I4,I4))+...
                    kron(kron(kron(M_mass(I1,I1),Stiff_1D(I2,I2)),M_mass(I3,I3)),M_mass(I4,I4))+...
                    kron(kron(kron(M_mass(I1,I1),M_mass(I2,I2)),Stiff_1D(I3,I3)),M_mass(I4,I4))+...
                    kron(kron(kron(M_mass(I1,I1),M_mass(I2,I2)),M_mass(I3,I3)),Stiff_1D(I4,I4));
                
                [II,JJ]=meshgrid(Real_index(Index_I+1),Real_index(Index_I+1));
                A_s=A_s+sparse(double(II),double(JJ),tmp,dof_sparse,dof_sparse);
                
                tmp=kron(kron(kron(b(I1),b(I2)),b(I3)),b(I4));
                b_s(Real_index(Index_I+1))=b_s(Real_index(Index_I+1))+tmp;
                uu_s(Real_index(Index_I+1))=uu_s(Real_index(Index_I+1))+...
                    kron(kron(kron(coef_MW(I1),coef_MW(I2)),coef_MW(I3)),coef_MW(I4));
                
                % save matrices to files
                A_encode{count}.A1=Stiff_1D(I1,I1);
                A_encode{count}.A2=M_mass(I2,I2);
                A_encode{count}.A3=M_mass(I3,I3);
                A_encode{count}.A4=M_mass(I4,I4);

                total_kron_flops = total_kron_flops + kron_mult_cost4( A1,A2,A3,A4);
                
                A_encode{count}.B1=M_mass(I1,I1);
                A_encode{count}.B2=Stiff_1D(I2,I2);
                A_encode{count}.B3=M_mass(I3,I3);
                A_encode{count}.B4=M_mass(I4,I4);

                total_kron_flops = total_kron_flops + kron_mult_cost4( B1,B2,B3,B4);
                
                A_encode{count}.C1=M_mass(I1,I1);
                A_encode{count}.C2=M_mass(I2,I2);
                A_encode{count}.C3=Stiff_1D(I3,I3);
                A_encode{count}.C4=M_mass(I4,I4);

                total_kron_flops = total_kron_flops + kron_mult_cost4(C1,C2,C3,C4);
                
                A_encode{count}.D1=M_mass(I1,I1);
                A_encode{count}.D2=M_mass(I2,I2);
                A_encode{count}.D3=M_mass(I3,I3);
                A_encode{count}.D4=Stiff_1D(I4,I4);

                total_kron_flops = total_kron_flops + kron_mult_cost4(D1,D2,D3,D4);
                
                A_encode{count}.IndexI=double(Real_index(Index_I+1));
                A_encode{count}.IndexJ=double(Real_index(Index_I+1));
                
                count=count+1;
                
                % Term 2: S*I*I*I--Assume j2_level==i2_level
                % j3_level==i3_level j4_level==i4_level
                j2_level=i2_level;j3_level=i3_level;j4_level=i4_level;
                for j1_level=0:i1_level-1
                    
                    J1=Index_1D(k,j1_level);
                    J2=I2;J3=I3;J4=I4;
                    
                    key_j=GenerateKey(J1,J2,J3,J4,Key1dMesh);
                    
                    Index_J=findvs(Hash,key_j);
                    
                    [II,JJ]=meshgrid(Real_index(Index_J+1),Real_index(Index_I+1));
                    
                    tmp=kron(kron(kron(Stiff_1D(J1,I1),M_mass(J2,I2)),M_mass(J3,I3)),M_mass(J4,I4));
                    
                    A_s=A_s+sparse([double(II),double(JJ)],[double(JJ),double(II)],[tmp',tmp'],dof_sparse,dof_sparse);
                    
                    % save matrices to files
                    A_encode{count}.A1=Stiff_1D(J1,I1);
                    A_encode{count}.A2=M_mass(J2,I2);
                    A_encode{count}.A3=M_mass(J3,I3);
                    A_encode{count}.A4=M_mass(J4,I4);

                    total_kron_flops = total_kron_flops + kron_mult_cost4(A1,A2,A3,A4);
                    
                    A_encode{count}.B1=0;
                    A_encode{count}.B2=0;
                    A_encode{count}.B3=0;
                    A_encode{count}.B4=0;
                    
                    A_encode{count}.C1=0;
                    A_encode{count}.C2=0;
                    A_encode{count}.C3=0;
                    A_encode{count}.C4=0;
                    
                    A_encode{count}.D1=0;
                    A_encode{count}.D2=0;
                    A_encode{count}.D3=0;
                    A_encode{count}.D4=0;
                    
                    A_encode{count}.IndexI=double(Real_index(Index_J+1));
                    A_encode{count}.IndexJ=double(Real_index(Index_I+1));
                    
                    A_encode{count+1}.A1=Stiff_1D(J1,I1)';
                    A_encode{count+1}.A2=M_mass(J2,I2)';
                    A_encode{count+1}.A3=M_mass(J3,I3)';
                    A_encode{count+1}.A4=M_mass(J4,I4)';

                    total_kron_flops = total_kron_flops + kron_mult_cost4(A1,A2,A3,A4);
                    
                    A_encode{count+1}.B1=0;
                    A_encode{count+1}.B2=0;
                    A_encode{count+1}.B3=0;
                    A_encode{count+1}.B4=0;
                    
                    A_encode{count+1}.C1=0;
                    A_encode{count+1}.C2=0;
                    A_encode{count+1}.C3=0;
                    A_encode{count+1}.C4=0;
                    
                    A_encode{count+1}.D1=0;
                    A_encode{count+1}.D2=0;
                    A_encode{count+1}.D3=0;
                    A_encode{count+1}.D4=0;
                    
                    A_encode{count+1}.IndexI=double(Real_index(Index_I+1));
                    A_encode{count+1}.IndexJ=double(Real_index(Index_J+1));
                    
                    count=count+2;
                    
                end
                
                % Term 3: I*S*I*I--Assume j1_level==i1_level
                % j3_level==i3_level j4_level==i4_level
                j1_level=i1_level;j3_level=i3_level;j4_level=i4_level;
                for j2_level=0:i2_level-1
                    
                    J2=Index_1D(k,j2_level);
                    J1=I1;J3=I3;J4=I4;
                    
                    key_j=GenerateKey(J1,J2,J3,J4,Key1dMesh);
                    
                    Index_J=findvs(Hash,key_j);
                    
                    [II,JJ]=meshgrid(Real_index(Index_J+1),Real_index(Index_I+1));
                    
                    tmp=kron(kron(kron(M_mass(J1,I1),Stiff_1D(J2,I2)),M_mass(J3,I3)),M_mass(J4,I4));
                    
                    A_s=A_s+sparse([double(II),double(JJ)],[double(JJ),double(II)],[tmp',tmp'],dof_sparse,dof_sparse);
                    
                    % save matrices to files
                    A_encode{count}.A1=0;
                    A_encode{count}.A2=0;
                    A_encode{count}.A3=0;
                    A_encode{count}.A4=0;
                    
                    A_encode{count}.B1=M_mass(J1,I1);
                    A_encode{count}.B2=Stiff_1D(J2,I2);
                    A_encode{count}.B3=M_mass(J3,I3);
                    A_encode{count}.B4=M_mass(J4,I4);

                    total_kron_flops = total_kron_flops + kron_mult_cost4(B1,B2,B3,B4);
                    
                    A_encode{count}.C1=0;
                    A_encode{count}.C2=0;
                    A_encode{count}.C3=0;
                    A_encode{count}.C4=0;
                    
                    A_encode{count}.D1=0;
                    A_encode{count}.D2=0;
                    A_encode{count}.D3=0;
                    A_encode{count}.D4=0;
                    
                    A_encode{count}.IndexI=double(Real_index(Index_J+1));
                    A_encode{count}.IndexJ=double(Real_index(Index_I+1));
                    
                    A_encode{count+1}.A1=0;
                    A_encode{count+1}.A2=0;
                    A_encode{count+1}.A3=0;
                    A_encode{count+1}.A4=0;
                    
                    A_encode{count+1}.B1=M_mass(J1,I1)';
                    A_encode{count+1}.B2=Stiff_1D(J2,I2)';
                    A_encode{count+1}.B3=M_mass(J3,I3)';
                    A_encode{count+1}.B4=M_mass(J4,I4)';

                    total_kron_flops = total_kron_flops + kron_mult_cost4(B1,B2,B3,B4);
                    
                    A_encode{count+1}.C1=0;
                    A_encode{count+1}.C2=0;
                    A_encode{count+1}.C3=0;
                    A_encode{count+1}.C4=0;
                    
                    A_encode{count+1}.D1=0;
                    A_encode{count+1}.D2=0;
                    A_encode{count+1}.D3=0;
                    A_encode{count+1}.D4=0;
                    
                    A_encode{count+1}.IndexI=double(Real_index(Index_I+1));
                    A_encode{count+1}.IndexJ=double(Real_index(Index_J+1));
                    count=count+2;
                    
                end
                
                % Term 4: I*I*S*I--Assume j1_level==i1_level
                % j2_level==i2_level j4_level==i4_level
                j1_level=i1_level;j2_level=i2_level;j4_level=i4_level;
                for j3_level=0:i3_level-1
                    
                    J3=Index_1D(k,j3_level);
                    J1=I1;J2=I2;J4=I4;
                    
                    key_j=GenerateKey(J1,J2,J3,J4,Key1dMesh);
                    
                    Index_J=findvs(Hash,key_j);
                    
                    [II,JJ]=meshgrid(Real_index(Index_J+1),Real_index(Index_I+1));
                    
                    tmp=kron(kron(kron(M_mass(J1,I1),M_mass(J2,I2)),Stiff_1D(J3,I3)),M_mass(J4,I4));
                    A_s=A_s+sparse([double(II),double(JJ)],[double(JJ),double(II)],[tmp',tmp'],dof_sparse,dof_sparse);
                    
                    % save matrices to files
                    A_encode{count}.A1=0;
                    A_encode{count}.A2=0;
                    A_encode{count}.A3=0;
                    A_encode{count}.A4=0;
                    
                    A_encode{count}.B1=0;
                    A_encode{count}.B2=0;
                    A_encode{count}.B3=0;
                    A_encode{count}.B4=0;
                    
                    A_encode{count}.C1=M_mass(J1,I1);
                    A_encode{count}.C2=M_mass(J2,I2);
                    A_encode{count}.C3=Stiff_1D(J3,I3);
                    A_encode{count}.C4=M_mass(J4,I4);

                    total_kron_flops = total_kron_flops + kron_mult_cost4(C1,C2,C3,C4);
                    
                    A_encode{count}.D1=0;
                    A_encode{count}.D2=0;
                    A_encode{count}.D3=0;
                    A_encode{count}.D4=0;
                    
                    A_encode{count}.IndexI=double(Real_index(Index_J+1));
                    A_encode{count}.IndexJ=double(Real_index(Index_I+1));
                    
                    A_encode{count+1}.A1=0;
                    A_encode{count+1}.A2=0;
                    A_encode{count+1}.A3=0;
                    A_encode{count+1}.A4=0;
                    
                    A_encode{count+1}.B1=0;
                    A_encode{count+1}.B2=0;
                    A_encode{count+1}.B3=0;
                    A_encode{count+1}.B4=0;
                    
                    A_encode{count+1}.C1=M_mass(J1,I1)';
                    A_encode{count+1}.C2=M_mass(J2,I2)';
                    A_encode{count+1}.C3=Stiff_1D(J3,I3)';
                    A_encode{count+1}.C4=M_mass(J4,I4)';

                    total_kron_flops = total_kron_flops + kron_mult_cost4(C1,C2,C3,C4);
                    
                    A_encode{count+1}.D1=0;
                    A_encode{count+1}.D2=0;
                    A_encode{count+1}.D3=0;
                    A_encode{count+1}.D4=0;
                    
                    A_encode{count+1}.IndexI=double(Real_index(Index_I+1));
                    A_encode{count+1}.IndexJ=double(Real_index(Index_J+1));
                    
                    count=count+2;
                end
                
                % Term 5: I*I*I*S--Assume j1_level==i1_level
                % j2_level==i2_level j3_level==i3_level
                j1_level=i1_level;j2_level=i2_level;j3_level=i3_level;
                for j4_level=0:i4_level-1
                    
                    J4=Index_1D(k,j4_level);
                    J1=I1;J2=I2;J3=I3;
                    
                    key_j=GenerateKey(J1,J2,J3,J4,Key1dMesh);
                    
                    Index_J=findvs(Hash,key_j);
                    
                    [II,JJ]=meshgrid(Real_index(Index_J+1),Real_index(Index_I+1));
                    
                    tmp=kron(kron(kron(M_mass(J1,I1),M_mass(J2,I2)),M_mass(J3,I3)),Stiff_1D(J4,I4));
                    A_s=A_s+sparse([double(II),double(JJ)],[double(JJ),double(II)],[tmp',tmp'],dof_sparse,dof_sparse);
                    
                    % save matrices to files
                    A_encode{count}.A1=0;
                    A_encode{count}.A2=0;
                    A_encode{count}.A3=0;
                    A_encode{count}.A4=0;
                    
                    A_encode{count}.B1=0;
                    A_encode{count}.B2=0;
                    A_encode{count}.B3=0;
                    A_encode{count}.B4=0;
                    
                    A_encode{count}.C1=0;
                    A_encode{count}.C2=0;
                    A_encode{count}.C3=0;
                    A_encode{count}.C4=0;
                    
                    A_encode{count}.D1=M_mass(J1,I1);
                    A_encode{count}.D2=M_mass(J2,I2);
                    A_encode{count}.D3=M_mass(J3,I3);
                    A_encode{count}.D4=Stiff_1D(J4,I4);

                    total_kron_flops = total_kron_flops + kron_mult_cost4(D1,D2,D3,D4);
                    
                    A_encode{count}.IndexI=double(Real_index(Index_J+1));
                    A_encode{count}.IndexJ=double(Real_index(Index_I+1));
                    
                    A_encode{count+1}.A1=0;
                    A_encode{count+1}.A2=0;
                    A_encode{count+1}.A3=0;
                    A_encode{count+1}.A4=0;
                    
                    A_encode{count+1}.B1=0;
                    A_encode{count+1}.B2=0;
                    A_encode{count+1}.B3=0;
                    A_encode{count+1}.B4=0;
                    
                    A_encode{count+1}.C1=0;
                    A_encode{count+1}.C2=0;
                    A_encode{count+1}.C3=0;
                    A_encode{count+1}.C4=0;
                    
                    A_encode{count+1}.D1=M_mass(J1,I1)';
                    A_encode{count+1}.D2=M_mass(J2,I2)';
                    A_encode{count+1}.D3=M_mass(J3,I3)';
                    A_encode{count+1}.D4=Stiff_1D(J4,I4)';

                    total_kron_flops = total_kron_flops + kron_mult_cost4(D1,D2,D3,D4);
                    
                    A_encode{count+1}.IndexI=double(Real_index(Index_I+1));
                    A_encode{count+1}.IndexJ=double(Real_index(Index_J+1));
                    
                    count=count+2;
                    
                end
                
                
            end
        end
    end
    
end
figure;
spy(A_s)

eigs(A_s,3,'SM')

disp(sprintf('nnz(A_s)=%g', nnz(A_s)));
disp(sprintf('total kron flops=%g ', total_kron_flops ));

% tic
sol_s = A_s\b_s*pi^2*4;
% toc

save(['./Data/A_4D_encode.mat'],'A_encode');


['Done of Solution']
% check error
norm(sol_s-uu_s)

% disp('Done of MapMW2DG')
%
% val=MapMW2DG'*sol_s;
% [x1,x2,x3,x4]=ndgrid(xxnode);
%
% sol=exactu_4D(x1(:),x2(:),x3(:),x4(:));
% norm(val-sol(:))
% max(abs(val-sol(:)))
% run ComputeError4D.m

function Ix=Index_1D(k,level)

if level==0
    Ix=[1:k];
else
    Ix=[k*2^(level-1)+1:k*2^level];
end

end

function key=GenerateKey(I1,I2,I3,I4,Key1dMesh)

tmp_1=Key1dMesh(I1,:);
tmp_2=Key1dMesh(I2,:);
tmp_3=Key1dMesh(I3,:);
tmp_4=Key1dMesh(I4,:);

n1=size(tmp_1,1);
n2=size(tmp_2,1);
n3=size(tmp_3,1);
n4=size(tmp_4,1);
dim=4;
ntol=n1*n2*n3*n4;
key=zeros(ntol,3*dim);
count=1;
for i1=1:n1
    for i2=1:n2
        for i3=1:n3
            for i4=1:n4
                key(count,[1:dim:3*dim])=tmp_1(i1,:);
                key(count,[2:dim:3*dim])=tmp_2(i2,:);
                key(count,[3:dim:3*dim])=tmp_3(i3,:);
                key(count,[4:dim:3*dim])=tmp_4(i4,:);
                count=count+1;
            end
        end
    end
end

key(:,3*dim-dim+1:3*dim)=key(:,3*dim-dim+1:3*dim)-1;
end

function key=GenerateKeyd(Id,Key1dMesh)

dim=size(Id,1);
n=zeros(dim,1);

ntol=1;
for i=1:dim
    tmp{:,i}=Key1dMesh(Id{i},:);
    n(i)=size(tmp{:,i},1);
    ntol=ntol*n(i);
end
ntol
% key=zeros(ntol,3*dim);
count=1;

for dd=dim:-1:1
    for i=1:n(dd)
        key(count,[dd:dim:3*dim])=tmp{:,dd}(i,:);
        count=count+1;
    end
end

% key(:,3*dim-dim+1:3*dim)=key(:,3*dim-dim+1:3*dim)-1;

end
%
% function [key,count]=recurseloop(Id,count)
% dim=size(Id,1);
% if count=
%     key(count,[dd])
%     count=count+1;
% end