function Y = kron_mult2(A1,A2,X)
% Y = kron_mult2(A1,A2,X)
%
Acell{1} = A1;
Acell{2} = A2;
nkron = 2;
Y = kron_multd(nkron,Acell,X);
