function pde=Vlasov4
% Numerical Example for Vlasov Equation
% Bump-on-tail instability


% parameters
k_0=0.3;A=0.04;Lmax=20*pi/3;Vmax=13;



    function f=Fx_0(x)
        % Initial condition for x variable
        f=(1+A*cos(k_0*x));
    end

    function f=Fv_0(v)
        % Initial condition for v variable
        np=9/(10*sqrt(2*pi));
        nb=2/(10*sqrt(2*pi));
        u=4.5;
        vt=0.5;
        f=np*exp(-v.^2/2)+nb*exp(-(v-u).^2/(2*vt^2));
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