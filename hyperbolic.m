function g=hyperbolic(h,tau,method_name)
gl=floor(1/h);
gw=floor(0.5/tau);
g=zeros(gl,gw);
lambda=tau/h;
tvd_vector=zeros(gw,1);

for i=1:gl
    g(i,1)=integral(@ini_f,h*(i-1),h*i)/h;
    if(i>1)
    tvd_vector(1)=tvd_vector(1)+abs((g(i,1)-g(i-1,1)));
    end
end
    function y=Riemann_Problem(uL,uR)
        if(abs(uL-uR)<1e-7)
            y=uL;
        elseif(uL>uR)
            y=uL;%sound speed >=0
        elseif(uL<uR && uL>=0)
            y=uL;
        elseif(uL<uR && uR<=0)
            y=uR;
        elseif(uL<=0 && uR>=0)
            y=0;
        end
    end
    function y=get_a(j)
        [q1,q2]=set_b1b2(j);
        y=(q1+q2)/2;
    end
    function [b1,b2]=set_b1b2(j)
        if(j==0)
            b1=0;
        else
            b1=g(j,n-1);
        end
        if(j==gl)
            b2=0;
        else
            b2=g(j+1,n-1);            
        end
    end
    function y=g_H(j)
        a=get_a(j);
        [m1,m2]=set_b1b2(j);
        y=(f(m1)+f(m2))/2-lambda*a*(f(m2)-f(m1))/2;
    end
    function y=g_L(j)
        [r1,r2]=set_b1b2(j);
        y=f(r1);
    end
    function y=g_flux(j)
        y=g_L(j)+phi(theta(j))*(g_H(j)-g_L(j));
    end
    function y=theta(j)
        a=get_a(j);
        [k1,k2]=set_b1b2(j);
        if(abs(k1-k2)<1e-7)
            y=0;
        elseif(a>=0)
            if(j==0)
                y=0;
            else
            [b3,b4]=set_b1b2(j-1);
            y=(b4-b3)/(k2-k1);
            end
        else
            if(j==gl)
                y=0;
            else
            [b3,b4]=set_b1b2(j+1);
            y=(b4-b3)/(k2-k1);
            end
        end
    end
f=@(u) u.^2/2;     
phi=@(theta) (abs(theta)+theta)./(1+abs(theta));
for n=2:gw
    if(string(method_name)==string('lax_wandroff'))
    for j=1:gl
        g(j,n)=g(j,n-1)-lambda*(g_H(j)-g_H(j-1));
    end
    elseif(string(method_name)==string('Godunov'))
    g(1,n)=g(1,n-1)-lambda*(f(Riemann_Problem(g(1,n-1),g(2,n-1)))-f(Riemann_Problem(0,g(1,n-1))));
    g(gl,n)=g(gl,n-1)-lambda*(f(Riemann_Problem(g(gl,n-1),0))-f(Riemann_Problem(g(gl-1,n-1),g(gl,n-1))));
    for j=2:(gl-1)
        g(j,n)=g(j,n-1)-lambda*(f(Riemann_Problem(g(j,n-1),g(j+1,n-1)))-f(Riemann_Problem(g(j-1,n-1),g(j,n-1))));
    end
    elseif(string(method_name)==string('flux_limiter'))
    for j=1:gl
        g(j,n)=g(j,n-1)-lambda*(g_flux(j)-g_flux(j-1));
    end
    elseif(string(method_name)==string('lax_friedrich'))
    for j=1:gl
        g(j,n)=g(j,n-1)-lambda*(g_L(j)-g_L(j-1));
    end
    elseif(string(method_name)==string('upwind'))
    for j=1:gl
                [b1,b2]=set_b1b2(j-1);
        g(j,n)=g(j,n-1)-lambda*(f(b2)-f(b1));
    end               

    end
        
end
end

