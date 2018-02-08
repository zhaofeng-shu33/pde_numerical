function [g,tvd_vector,mse]=hyperbolic(h,tau,method_name)
L=floor(2/h);
grid_width=floor(1/tau);
g=zeros(L,grid_width);
lambda=tau/h;
tvd_vector=zeros(grid_width,1);

for i=1:L
    g(i,1)=integral(@ini_f,h*(i-1),h*i)/h;
    if(i>1)
    tvd_vector(1)=tvd_vector(1)+abs((g(i,1)-g(i-1,1)));
    end
end

for n=2:grid_width
    if(string(method_name)==string('lax_friedrich'))
    g(1,n)=0.5*(1-lambda)*g(2,n-1);
    g(L,n)=0.5*(1+lambda)*g(L-1,n-1);
    for j=2:(L-1)
        g(j,n)=0.5*(g(j+1,n-1)+g(j-1,n-1))-0.5*lambda*(g(j+1,n-1)-g(j-1,n-1));
        tvd_vector(n)=tvd_vector(n)+abs((g(j,n)-g(j-1,n)));
    end
      tvd_vector(n)=tvd_vector(n)+abs((g(L,n)-g(L-1,n)));
    elseif(string(method_name)==string('upwind'))
    for j=2:(L)
        g(j,n)=g(j,n-1)-lambda*(g(j,n-1)-g(j-1,n-1));
        tvd_vector(n)=tvd_vector(n)+abs((g(j,n)-g(j-1,n)));
    end
    elseif(string(method_name)==string('lax_wandroff'))
    g(1,n)=g(1,n-1)-0.5*lambda*g(2,n-1)+0.5*lambda^2*(g(2,n-1)-2*g(1,n-1));
    g(L,n)=g(L,n-1)+0.5*lambda*g(L-1,n-1)+0.5*lambda^2*(g(L-1,n-1)-2*g(L,n-1));
    for j=2:(L-1)
        g(j,n)=g(j,n-1)-lambda/2*(g(j+1,n-1)-g(j-1,n-1))+lambda^2/2*(g(j+1,n-1)-2*g(j,n-1)+g(j-1,n-1));
        tvd_vector(n)=tvd_vector(n)+abs((g(j,n)-g(j-1,n)));
    end
        tvd_vector(n)=tvd_vector(n)+abs((g(L,n)-g(L-1,n)));
    elseif(string(method_name)==string('beam_warming'))
    g(1,n)=(1-1.5*lambda+lambda^2/2)*g(1,n-1);
    g(2,n)=(1-1.5*lambda+lambda^2/2)*g(2,n-1)+(2*lambda-lambda^2)*g(1,n-1);
    tvd_vector(n)=tvd_vector(n)+abs((g(2,n))-g(1,n));
    for j=3:L
        g(j,n)=g(j,n-1)-0.5*lambda*(3*g(j,n-1)-4*g(j-1,n-1)+g(j-2,n-1))+0.5*lambda^2*(g(j,n-1)-2*g(j-1,n-1)+g(j-2,n-1));
        tvd_vector(n)=tvd_vector(n)+abs((g(j,n)-g(j-1,n)));
    end
       
    end
        

end
x=linspace(0,2-h,2/h)+h/2;
y=ini_f(x-1);
mse=sum(power(y-g(:,end)',2));
end

