function isok = fk6d_vlasov7_check(testCase)
% isok = fk6d_vlasov7_check(testCase)
%
% Test for PDE case 7, which is coupled with given E
% compression
% 0:: Explicitly construct the A matrix and store all non-zeros, apply as A*f, no tensor product encoding
% 1:: Explicitly construct the A matrix in sparse storage format, apply as A*f, no tensor product encoding
% 2:: Do not construct A, only create elements as needed from the A_data data structure, apply as the element-wise breakdown of A*f, no tensor product encoding (compression=2).
% 3:: Do not construct A, use A_encode, apply with kron_multd, tensor product encoding of the Deg DOFs
% 4:: Do not construct A, use A_data, apply with kron_multd, tensor product encoding of the Deg DOFs

  idebug = 1;
  addpath(genpath(pwd));
  disp(sprintf('Testing with valsov10, testCase=%d',testCase));
  
  switch(testCase)
   case 1 
      quiet = 1; lev = 3; deg = 2; TEND = 1; compression = 0;
   case 2 
      quiet = 1; lev = 3; deg = 2; TEND = 1; compression = 1;
   case 3 
      quiet = 1; lev = 3; deg = 2; TEND = 1; compression = 2;
   case 4 
      quiet = 1; lev = 3; deg = 2; TEND = 1; compression = 3;
   case 5 
      quiet = 1; lev = 3; deg = 2; TEND = 1; compression = 4;
   otherwise 
      error(sprintf('fk6d_vlasov7_check: invalid testCase=%d',testCase));
  end;
  
  
  act_f = fk6d(Vlasov7,lev,deg,TEND,quiet,compression);
  

  load('tests/vlasov7/solution.mat');
  exp_f = fval;
  
  isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
  
  tol = 1e-4;
  isok = all(   abs(act_f(:)-exp_f(:)) <= tol * abs( exp_f(:) ) );

  if (idebug >= 1),
    if (~isok),
       maxerr = max(abs(act_f(:)-exp_f(:)));
       error(sprintf('testCase=%d, maxerr = %g', testCase,maxerr ));
    end;
  end;


  if (~isOctave),
    verifyEqual(testCase,act_f,exp_f,'RelTol',tol);
  end;

end 

