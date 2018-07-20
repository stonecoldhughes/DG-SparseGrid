%alternating flux
A1=[2.0000e+00   3.0000e+00   1.0000e-02   3.8491e-01   1.0291e+00
    2.0000e+00   4.0000e+00   7.0711e-03   1.0620e-01   3.0045e-01
    2.0000e+00   5.0000e+00   5.0000e-03   4.0002e-02   1.5498e-01
    2.0000e+00   6.0000e+00   3.5355e-03   5.6876e-03   6.0936e-02
    2.0000e+00   7.0000e+00   2.5000e-03   2.3206e-03   2.1275e-02
    2.0000e+00   8.0000e+00   1.7678e-03   5.9037e+32   7.5693e+33];%cfl=0.02

A2=[2.0000e+00   3.0000e+00   1.5000e-02   3.8209e-01   1.0203e+00
    2.0000e+00   4.0000e+00   1.0607e-02   1.0526e-01   2.9770e-01
    2.0000e+00   5.0000e+00   7.5000e-03   3.9869e-02   1.5435e-01
    2.0000e+00   6.0000e+00   5.3033e-03   5.7724e-03   6.1128e-02
    2.0000e+00   7.0000e+00   3.7500e-03   5.1447e+26   4.7126e+27];%cfl=0.03

A3=[2.0000e+00   3.0000e+00   5.0000e-03   3.8459e-01   1.0301e+00
    2.0000e+00   4.0000e+00   3.5355e-03   1.0782e-01   3.0074e-01
    2.0000e+00   5.0000e+00   2.5000e-03   4.0311e-02   1.5642e-01
    2.0000e+00   6.0000e+00   1.7678e-03   5.7984e-03   6.1219e-02
    2.0000e+00   7.0000e+00   1.2500e-03   1.5378e-03   2.0851e-02
    2.0000e+00   8.0000e+00   8.8388e-04   8.4473e-04   7.1600e-03];%cfl=0.01

for i =2:5
    a1(i,1)=log2(A1(i-1,3)/A(i,3));
    a2(i,1)=log2(A1(i-1,4)/A(i,4));
    c1(i,1)=log2(A2(i-1,3)/C(i,3));
    c2(i,1)=log2(A2(i-1,4)/C(i,4));
    u1(i,1)=log2(A3(i-1,3)/U(i,3));
    u2(i,1)=log2(A3(i-1,4)/U(i,4));
end

