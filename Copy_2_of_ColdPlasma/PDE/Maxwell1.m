function pde=Maxwell1
%====================================
% Test 1 for time harmonic
% Maxwells' Equations
%====================================


mu = 1;
eps = 1;
w2 = 1e-0;

    function E=exact_E(x,y,z)
        E(:,1)= sin(pi*y).*sin(pi*z);
        E(:,2)= sin(pi*z).*sin(pi*x);
        E(:,3)= sin(pi*x).*sin(pi*y);
    end


    function f=rhs(x,y,z)
        f(:,1) = (2*pi^2-pde.w2)*sin(pi*y).*sin(pi*z);
        f(:,2) = (2*pi^2-pde.w2)*sin(pi*z).*sin(pi*x);
        f(:,3) = (2*pi^2-pde.w2)*sin(pi*x).*sin(pi*y);
        
        f = f/(2*pi^2-pde.w2);
    end



pde = struct('mu',mu,'eps',eps,'w2',w2,'E',@exact_E,'rhs',@rhs);

end