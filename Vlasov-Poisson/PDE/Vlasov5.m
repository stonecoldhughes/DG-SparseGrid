function pde=Vlasov5
% Numerical Example for Vlasov Equation
% Two-stream instability

% parameters
k_0=0.5;A=0.05;Lmax=4*pi;Vmax=2*pi;



    function f=Fx_0(x)
        % Initial condition for x variable
        f=(1+A*cos(k_0*x));
    end

    function f=Fv_0(v)
        % Initial condition for v variable
        f=v.^2.*exp(-v.^2/2)/sqrt(2*pi);
    end
    function f=Fxv_0(x,v)
        f=Fv_0(v).*Fx_0(x);
    end
    function f=Ex(x)
        f=x-x;
    end
    function f=Et(x)
        f=x-x;
        
    end
    function f=E(x,t)
%         f=x-x;
%         f=A/k_0*exp(-k_0^2*t^2/2).*sin(k_0*x)-x+8*pi;
            f=x-x;%A/k_0*sin(k_0*x).*exp(-k_0^2*t.^2/2)-x+8*pi;
    end
    function f=F(x,v,t)
        f=Fv_0(v).*(A*cos(k_0*(x-v*t)));
    end
    function f=rho(x,t)
%        f=exp(-k_0*t.^2/2).*cos(k_0*x); 
        f= -A*cos(k_0*x).*exp(-k_0*t);
    end
    
pde = struct('k_0',k_0,'Fx_0',@Fx_0,'Fv_0',@Fv_0,'Fxv_0',@Fxv_0,'Ex',@Ex,...
    'Et',@Et,'E',@E,'F',@F,'rho',@rho,'Vmax',Vmax,'Lmax',Lmax);

end