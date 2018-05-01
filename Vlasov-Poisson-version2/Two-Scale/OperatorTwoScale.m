function FMWT_COMP = OperatorTwoScale(maxDeg,maxLev)
%----------------------------------
% Set-up Two-scale operator       %
%----------------------------------
% Input: Degree: maxDeg
%        Level: Np
% Output: Convert Matrix: FMWT_COMP
%**********************************

% Load G0 and H0 from file

load(['two_scale_rel_',num2str(maxDeg),'.mat'])


H0(find(abs(H0)<1e-5))=0; % Why are we doing this?
G0(find(abs(G0)<1e-5))=0;

H1 = zeros(maxDeg);
G1 = zeros(maxDeg);

for j_x = 1:maxDeg
    for j_y = 1:maxDeg
        H1(j_x,j_y) = ((-1)^(j_x+j_y-2)  )*H0(j_x,j_y);
        G1(j_x,j_y) = ((-1)^(maxDeg+j_x+j_y-2))*G0(j_x,j_y);
    end
end

FMWT = zeros(maxDeg*maxLev);

for j=1:maxLev/2
    % The reverse order from Lin
    FMWT(maxDeg*(j-1)+1:maxDeg*j,2*maxDeg*(j-1)+1:2*maxDeg*j)=[H0 H1];
    FMWT(maxDeg*(j+maxLev/2-1)+1:maxDeg*(j+maxLev/2),2*maxDeg*(j-1)+1:2*maxDeg*j) = [G0 G1];
end

FMWT_COMP = eye(maxDeg*maxLev);

n=log2(maxLev);

for j=1:n
    cFMWT = FMWT;
    % Reverse the index in matrix from Lin
    if j>1
        cFMWT = zeros(maxDeg*maxLev);
        cn = 2^(n-j+1)*maxDeg;
        cnr=maxLev*maxDeg-cn;
        cFMWT(cn+1:maxDeg*maxLev,cn+1:maxDeg*maxLev)=eye(maxLev*maxDeg-cn);
        cFMWT(1:cn/2,1:cn)=FMWT(1:cn/2,1:cn);
        cFMWT(cn/2+1:cn,1:cn)=FMWT(maxDeg*maxLev/2+1:maxDeg*maxLev/2+cn/2,1:cn);
    end

    FMWT_COMP = cFMWT*FMWT_COMP;
end

