%Modelo Dengue - Fracionário
dt=0.001;		% espacamento de tempo
tf =365;	    % tempo final
t=0:dt:tf;		% tempo	0 ate tempo final
t0=0;
% Condicoes iniciais
temp=(11.19/2)*cos((2*pi.*t)/365)+22.77; %Função de temperatura ao longo do ano
plot(t,temp)                            
prec=(63.2/2)*cos((2*pi.*t)/365)+(63.2/2);%Função de precipitação ao longo do ano
HS_0=379138;
HI_0=8;
HR_0=0;
MS_0=568719;
MI_0=0;
A_0=0;
muh=3.9E-5;
sigma=(1/7);
H=379146; %População Bauru
k=0.8;
delta=-15.837+1.2897.*temp-0.0163.*temp.^2;
plot(t,delta)
b=0.056.*delta;   
plot(t,b)
gamma=1;
tau=1;
betah=0.4;
betam=0.033.*temp-0.41;
plot(t,betam)
mum=0.8962-0.159.*temp+(1.116E-2).*temp.^2-(3.408E-4).*temp.^3+(3.809E-6).*temp.^4;
plot(t,mum)
mua=(2.13-0.3797.*temp+(2.457E-2).*temp.^2-(6.778E-4).*temp.^3+(6.794E-6).*temp.^4)/7;
plot(t,mua)
alpha=(0.131-(5.723E-2).*temp+(1.164E-2).*temp.^2-(1.341E-3).*temp.^3+(8.723E-5).*temp.^4-(3.017E-6).*temp.^5+(5.153E-8).*temp.^6+(3.42E-10).*temp.^7)/7;
plot(t,alpha)
Cmax=1000;
C=(Cmax/63.2).*prec;
plot(t,C)

f=@(t,y)[(tau.^(1-gamma))*(muh*(H-y(1))-(b*betah*y(1)*y(5))/(H));
    (tau.^(1-gamma))*((b*betah*y(1)*y(5))/(H)-(muh+sigma)*y(2));
    (tau.^(1-gamma))*(sigma*y(2)-muh*y(3));
    (tau.^(1-gamma))*(alpha*y(6)-(b*betam*y(4)*y(2))/(H)-mum*y(4));
    (tau.^(1-gamma))*((b*betam*y(4)*y(2))/(H)-mum*y(5));
    (tau.^(1-gamma))*((k*delta*(1-(y(6)/C))*(y(4)+y(5)))-(mua+alpha)*y(6))]; %

[f,y]=fde12(gamma,f,t0,tf,[HS_0;HI_0;HR_0;MS_0;MI_0;A_0],dt);			

HS_estimado=y(1,:);				% Soluções
HI_estimado=y(2,:);	
HR_estimado=y(3,:);	
MS_estimado=y(4,:);	
MI_estimado=y(5,:);	
A_estimado=y(6,:);	