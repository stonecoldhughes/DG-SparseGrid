%Fixed dt
%dt=1/1000,MaxT=100
C0=[2.0000e+00   4.0000e+00   1.2985e-02   3.9679e-02
    2.0000e+00   5.0000e+00   4.9623e-03   2.0154e-02
    2.0000e+00   6.0000e+00   1.7881e-03   1.0117e-02
    2.0000e+00   7.0000e+00   6.3519e-04   5.0646e-03
    2.0000e+00   8.0000e+00   2.2484e-04   2.5353e-03
    2.0000e+00   9.0000e+00   7.9516e-05   1.2725e-03
    2.0000e+00   1.0000e+01   7.5191e-05   6.4563e-04
    2.0000e+00   1.1000e+01   7.4913e-05   3.4096e-04];

%dt=1/100,MaxT=10
C1=[2.0000e+00   4.0000e+00   1.2087e-02   3.7288e-02
    2.0000e+00   5.0000e+00   4.6132e-03   1.8940e-02
    2.0000e+00   6.0000e+00   1.6618e-03   9.5656e-03
    2.0000e+00   7.0000e+00   7.4749e-04   4.9101e-03
    2.0000e+00   8.0000e+00   7.2935e-04   2.6879e-03
    2.0000e+00   9.0000e+00   7.2482e-04   1.7324e-03
    2.0000e+00   1.0000e+01   7.2369e-04   1.3947e-03
    2.0000e+00   1.1000e+01   7.2340e-04   1.2966e-03];

%Fixed CFL Time period=1
%CFL=1                        %condition number
C2=[2.0000e+00   4.0000e+00   4.0932e+01   8.2119e-02   1.4951e-01
    2.0000e+00   5.0000e+00   6.6685e+01   4.1434e-02   7.5300e-02
    2.0000e+00   6.0000e+00   7.8923e+01   2.0284e-02   3.7700e-02
    2.0000e+00   7.0000e+00   7.5590e+01   9.9590e-03   1.8888e-02
    2.0000e+00   8.0000e+00   8.7569e+01   4.9297e-03   9.4899e-03
    2.0000e+00   9.0000e+00   9.6451e+01   2.4526e-03   4.7724e-03
    2.0000e+00   1.0000e+01   8.9428e+01   1.2233e-03   2.3968e-03
    2.0000e+00   1.1000e+01   9.3258e+01   6.1092e-04   1.2016e-03];
%Lev=11 cost 414s, inv(H) 100s

%CFL=0.1                      %condition number
C3=[2.0000e+00   4.0000e+00   6.1163e+00   4.0492e-02   1.0394e-01
    2.0000e+00   5.0000e+00   7.2906e+00   1.5574e-02   5.2493e-02
    2.0000e+00   6.0000e+00   7.4859e+00   5.6453e-03   2.6501e-02
    2.0000e+00   7.0000e+00   6.9434e+00   2.0110e-03   1.3336e-02
    2.0000e+00   8.0000e+00   8.1802e+00   7.1249e-04   6.6919e-03
    2.0000e+00   9.0000e+00   7.6361e+00   2.5204e-04   3.3522e-03
    2.0000e+00   1.0000e+01   8.5874e+00   1.2334e-04   1.6777e-03
    2.0000e+00   1.1000e+01   8.8967e+00   6.1347e-05   8.3926e-04];

%CFL=3
C4=[2.0000e+00   4.0000e+00   8.7641e+01   2.2586e-01   3.0940e-01
    2.0000e+00   5.0000e+00   1.0893e+02   1.1807e-01   1.6484e-01
    2.0000e+00   6.0000e+00   1.3834e+02   6.1363e-02   8.4174e-02
    2.0000e+00   7.0000e+00   2.1189e+02   3.1148e-02   4.2867e-02
    2.0000e+00   8.0000e+00   2.3159e+02   1.4783e-02   2.1219e-02
    2.0000e+00   9.0000e+00   2.6963e+02   7.4691e-03   1.0797e-02
    2.0000e+00   1.0000e+01   1.7290e+02   3.6606e-03   5.4062e-03
    2.0000e+00   1.1000e+01   2.8076e+02   1.8398e-03   2.7241e-03];

%Lev=11 time cose for inv(H) is 80s

%plot numerical solution and exact solution 
%can not use period time=1, since cos(pi/2)=0



    

   


    