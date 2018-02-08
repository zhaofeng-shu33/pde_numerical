#Implicit RK symplectic Method
T<-matrix(,4,4)
h=0.5
J2=matrix(c(0,-1,1,0),2,2)
T[1:2,1:2]=diag(2)+h*J2/4
T[1:2,3:4]=(0.25-sqrt(3)/6)*h*J2
T[3:4,1:2]=(0.25+sqrt(3)/6)*h*J2
T[3:4,3:4]=diag(2)+h*J2/4
z=matrix(,2,as.integer(100/h))
z[,1]=c(0,1)
for(i in 2:length(z[1,])){
b=-c(J2%*%z[,i-1],J2%*%z[,i-1])
result=solve(T,b)
K1=result[1:2]
K2=result[3:4]
z[,i]=z[,i-1]+h*(K1+K2)/2
}
z_energy=(z[1,]^2+z[2,]^2)/2
#explicit,classical Runge-Kutta
y=matrix(,2,as.integer(100/h))
y[,1]=c(0,1)
for(i in 2:length(y[1,])){
K1=-J2%*%y[,i-1]
K2=-J2%*%(y[,i-1]+0.5*h*K1)
K3=-J2%*%(y[,i-1]+0.5*h*K2)
K4=-J2%*%(y[,i-1]+h*K3)
y[,i]=y[,i-1]+h*(K1/6+K2/3+K3/3+K4/6)
}
y_energy=(y[1,]^2+y[2,]^2)/2
