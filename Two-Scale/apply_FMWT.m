function Y = apply_FMWT( kdeg, Lev, FMWT, X, imethod_in )
%
% Y = apply_FMWT( k, lev, FMWT, X )
% 
% Y = FMWT * X
%
idebug = 0;
if (idebug >= 1),
        t1 = cputime();
end;
n = kdeg * 2^Lev;
isok  = (size(FMWT,1) == n) && (size(FMWT,2) == n) && ...
        (size(X,1) == n );
nvec = size(X,2);
if (~isok),
    error(sprintf('apply_FMWT: invalid sizes, size(FMWT)=(%d,%d),size(X)=(%d,%d)', ...
              size(FMWT,1), size(FMWT,2),    size(X,1),size(X,2)) );
    return;
end;

imethod_default = 3;

imethod = imethod_default;
if (nargin >= 5),
        imethod = imethod_in;
end;

if (imethod == 1),
      % ----------------------------
      % just perform matrix multiply
      % ----------------------------
      Y = FMWT * X;
elseif (imethod == 2),
      % -------------------------------------------- 
      % take advantage of special sparsity structure
      % -------------------------------------------- 
      Y = zeros( size(FMWT,1), size(X,2) );

      % -------------------------------
      % special case for coarsest level
      % -------------------------------
      ip = 1;
      ipend = 2*kdeg;
      col1 = 1;
      col2 = n;
      Y(ip:ipend,1:nvec) = FMWT(ip:ipend, col1:col2) * X( col1:col2,1:nvec);

      ip = 2*kdeg + 1;
      for lev=1:(Lev-1),
          ncells = 2^lev;
          isize = n/ncells;
          for icell=1:ncells,
             ipend = ip + kdeg-1;
                  col1 = 1 + (icell-1)*isize;
                  col2 = col1 + isize-1;
                  Y(ip:ipend,1:nvec) = FMWT(ip:ipend,  col1:col2) * X( col1:col2,1:nvec);
             ip = ipend + 1;
          end;
      end;
elseif (imethod == 3),

      % -------------------------------------------- 
      % take advantage of special sparsity structure
      % and identical sub-matrices
      % -------------------------------------------- 
      Y = zeros( size(FMWT,1), size(X,2) );

      % -------------------------------
      % special case for coarsest level
      % -------------------------------
      ip = 1;
      ipend = 2*kdeg;
      col1 = 1;
      col2 = n;
      Y(ip:ipend,1:nvec) = FMWT(ip:ipend, col1:col2) * X( col1:col2,1:nvec);

      ip = 2*kdeg + 1;
      for lev=1:(Lev-1),
          ncells = 2^lev;
          isize = n/ncells;

          % --------------------------------------
          % take advantage of identical sub-blocks
          % to make fewer calls to batched GEMM
          % --------------------------------------
          icell = 1;
          ipend = ip + kdeg-1;
          col1 = 1 + (icell-1)*isize;
          col2 = col1 + isize-1;

          Fmat = FMWT(ip:ipend,col1:col2);
          ip2 = ip + (ncells * isize)-1;

          ncol = (n*nvec/isize);
          nrows = (ncells*kdeg);
          ipend = ip + (ncells*kdeg)-1;

          Y(ip:ipend,1:nvec) = ...
           reshape(  Fmat(1:kdeg,1:isize)*reshape(X(1:n,1:nvec), isize, ncol), ...
                     nrows, (kdeg*ncol)/nrows  );


          ip = ip + (ncells * kdeg);

      end;


else
    error(sprintf('apply_FMWT: invalid imethod=%d',imethod));
    return;
end;

if (idebug >= 1),
   t2 = cputime();
   disp(sprintf('apply_FMWT: n=%d, nvec=%d, imethod=%d, time=%g', ...
        n,  nvec, imethod, (t2-t1) ));
end;

end

   

