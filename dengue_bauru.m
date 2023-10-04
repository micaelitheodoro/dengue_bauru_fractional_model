%Modelo Dengue - Fracionário
dt=0.001;		% espacamento de tempo
tf =365;	    % tempo final
t=0:dt:tf;		% tempo	0 ate tempo final
t0=0;
% Condicoes iniciais

HS_0=379138;
HI_0=8;
HR_0=0;
MS_0=568719;
MI_0=0;
A_0=0;
muh=3.9E-5;
b=1;
sigma=(1/7);
H=379146; %População Bauru
C=1000;
k=0.8;
delta=6.353;
%mua=0.061;       
gamma=1;
tau=1;
betah=0.4;
betam=0.4;
mum=x(1);
mua=x(2);
alpha=x(3);

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