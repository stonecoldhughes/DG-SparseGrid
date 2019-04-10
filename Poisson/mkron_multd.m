function Y = mkron_multd( nkron, Acell, X, offset_in )
%  Y = mkron_multd( nkron, Acell, X )
%
%  generic multi-dimension kron product
%
%  implementation to minimize flops
%
%  note using Acell{offset+1}, ... Acell{offset+nkron}
%
%
offset = 0;
if (nargin >= 4),
  offset = offset_in;
end;
n = nkron;

isok = (length(Acell) >= offset+nkron);
if (~isok),
  error(sprintf('kron_multd: invalid nkron=%d, length(Acell)=%d', ...
        nkron, length(Acell) ));
  return;
end;

  rc = zeros(2,nkron);
  for k=1:nkron,
    rc(1,k) = size( Acell{offset+k},1);
    rc(2,k) = size( Acell{offset+k},2);
  end;
  
  isizeX = prod( rc(2,1:nkron) );
  isizeY = prod( rc(1,1:nkron) );
  
  nvec = prod(size(X))/isizeX;
  isok = (mod( prod(size(X)), isizeX) == 0);
  if (~isok),
    error(sprintf('kron_multd:nkron=%d,prod(size(X))=%g,isizeX=%g', ...
               nkron, prod(size(X)), isizeX ));
  end;

if (nkron == 1),
  nc = size(Acell{offset+1},2);
  Y = Acell{offset+1} * reshape(X, [nc, prod(size(X))/nc]);
elseif (nkron == 2),
  [flops,isplit,imethod] = kron_minflops(rc(1:2,1:nkron));
  A = Acell{offset+1}; nrowA = size(A,1); ncolA = size(A,2);
  B = Acell{offset+2}; nrowB = size(B,1); ncolB = size(B,2);
  nvec = prod(size(X))/(ncolB*ncolA);
  X = reshape( X, ncolB*ncolA, nvec);
  Y = zeros( nrowB*nrowA, nvec);

  if (imethod == 1),
    % ----------------------
    % Y = (B*X)*transpose(A)
    % BX = B * X;
    % ----------------------
    
    X = reshape(X, ncolB, (ncolA*nvec));
    BX = reshape(B * X, [nrowB*ncolA,nvec]);
    
    for i=1:nvec,
       BXi = reshape( BX(:,i), nrowB, ncolA);
       Yi = BXi(1:nrowB,1:ncolA) * transpose( A(1:nrowA,1:ncolA)); % Yi is nrowB by nrowA
       Y(:,i) = reshape(Yi, nrowB*nrowA,1);
    end;
  else
    % ------------------------
    % Y = B * (X*transpose(A))
    % ------------------------
    XAt = zeros(ncolB*nrowA,nvec);
    for i=1:nvec,
      Xi = reshape(X(:,i), ncolB,ncolA);
      XiAt = Xi(1:ncolB,1:ncolA) * transpose(A(1:nrowA,1:ncolA)); % XiAt is ncolB by nrowA
      XAt(:,i) = reshape( XiAt, ncolB*nrowA,1);
    end;
    XAt = reshape( XAt, ncolB, nrowA*nvec );
    Y = B(1:nrowB,1:ncolB) * XAt(1:ncolB, 1:(nrowA*nvec) );
    Y = reshape( Y, nrowB*nrowA, nvec );
  end; % if (imethod)
     
elseif (nkron >= 3),     
  % ------------------------------
  % general case require recursion
  % ------------------------------
  [flops,isplit,imethod] = kron_minflops( rc ); 
  k = isplit;
  isok = (1 <= k) && (k <= (n-1));
  if (~isok),
    disp(sprintf('nkron=%d,k=%d ',nkron,k));
    disp(sprintf('flops=%g,isplit=%g,imethod=%g', ...
                  flops,   isplit,   imethod ));
    for i=1:size(rc,2),
       disp(sprintf('rc(1:2,%d)=(%d,%d)', ...
                            i, rc(1,i), rc(2,i) ));
    end;
  end;

  kp1 = k+1;
  n1 = prod(rc(2,1:k));
  n2 = prod(rc(2,kp1:n));
  m1 = prod(rc(1,1:k));
  m2 = prod(rc(1,kp1:n));

  X = reshape( X, n1*n2, nvec);
  Y = zeros( m1*m2,nvec );
  if (imethod == 1),
    %  X = reshape(X, n2,n1)
    %  T1 = kron( Akp1...An) *  X
    %  T1 is m2 by n1
    %  T1t = reshape( transpose(T1), n1, m2 )
    %  Y = transpose(kron( A1..Ak) * T1t)
  
    

       X = reshape( X, n2, n1*nvec );
       T1all = mkron_multd( (n-kp1+1), Acell, X, offset+k ); % note T1all is m2 by (n1*nvec)
       T1all = reshape( T1all, m2*n1, nvec );

      T1tall = zeros(n1*m2, nvec );
      for i=1:nvec,
       T1 = reshape( T1all(:,i), m2,n1);
       T1t = reshape( transpose(T1), n1,m2);
       T1tall(:,i) = reshape( T1t, n1*m2,1);
      end;

       Yitall = mkron_multd( k, Acell, T1tall,offset+0); % Y1tall is m1 by m2*nvec
       Yitall = reshape( Yitall, m1*m2, nvec);
       for i=1:nvec,
         Yit = reshape( Yitall(:,i), m1,m2);
         Yi =  reshape(transpose(Yit), m2,m1); 
         Y(1:(m1*m2),i) = reshape(Yi, m1*m2,1);
        end;
  else
    % Xt = reshape( transpose(X), n1,n2);
    % T2 = kron( A1..Ak ) * Xt, T2 is m1 by n2
    % T2t = reshape(transpose(T2), n2, m1 )
    % Y = kron( Akp1..An) * T2t
   
    Xt = zeros( n1*n2, nvec );
    for ivec=1:nvec,
      Xi = reshape(X(:,ivec), n2,n1);
      Xit = reshape(transpose(Xi), n1, n2);
      Xt(:,ivec) = reshape(Xit, n1*n2,1);
    end;
    Xt = reshape( Xt, n1, n2*nvec );
    T2all = mkron_multd( k, Acell, Xt,offset+0 ); % note T2all is m1 by (n2*nvec)
    isok = (prod(size(T2all)) == (m1*n2*nvec));
    if (~isok),
      disp(sprintf('size(T2all)=%d ~= m1*n2*nvec=%d', ...
                    prod(size(T2all)), m1*n2*nvec));
      disp(sprintf('m1=%d,m2=%d,n1=%d,n2=%d,nvec=%d,size(T2all)=%g',...
                    m1,   m2,   n1,   n2,   nvec,   prod(size(T2all)) ));
      disp(sprintf('k=%d,offset=%d',k,offset));
      for i=1:nkron,
        disp(sprintf('rc(1:2,%d)=(%d,%d)', ...
                             i,  rc(1,i), rc(2,i) ));
      end;
    end;

    T2all = reshape( T2all, m1*n2, nvec );

    T2tall = zeros( m1*n2, nvec );
    for i=1:nvec,
      T2 = reshape( T2all(:,i), m1, n2 );
      T2t = reshape( transpose(T2), n2, m1 );
      T2tall(:,i) = reshape( T2t, n2*m1,1);
    end;
    T2tall = reshape( T2tall, n2, m1*nvec );
      
      Yi = mkron_multd( (n-kp1+1), Acell, T2tall, offset+k ); % Yi is m2 by m1*nvec

      Y(1:(m1*m2),1:nvec) = reshape( Yi, m1*m2,nvec);
  end; % if (imethod)
 end; % if (nkron == 1)
end;