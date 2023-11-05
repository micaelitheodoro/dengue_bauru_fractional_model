%Modelo Dengue - Fracionário
%clc
%clear all
dt=0.0001;		
tf =365;	    
t=0:dt:tf;		
t0=0;
temp=((11.19/2)*cos((2*pi*t)/365)+22.77)'; %Função de temperatura ao longo do ano  
%plot(t,temp);
prec=((63.2/2)*cos((2*pi*t)/365)+(63.2/2))';%Função de precipitação ao longo do ano
HS_0=379138;
HI_0=8;
HR_0=0;
MS_0=568719;
MI_0=0;
A_0=0;
%muh=3.9E-5;
muh=0.0000346417; %taxa de mortalidade Bauru
sigma=(1/7);
H=379146;%População Bauru
% H = ones(3650001,1)*379146;
k=0.8;
delta=-15.837+(1.2897*temp)-0.0163*(temp).^2;%Taghikhani2018
%delta=-5.4+18*temp-0.2124*temp.^2+(1.015E-2)*temp.^3-(1.515E-4)*temp.^4;%Putra 2017
b=(0.056*delta);   %Lourdes 2015
%b=0.0943+0.0043*temp;%Taghikhani2018
% b=(0.03*temp+0.66)/7; %Putra
gamma=x(1);
tau=x(2);
% betah=0.4;
% betah=0.133707252279892;
%betah=x(1);
 betah=0.023*temp+0.122;%Lourdes 2015
%  for k=1:1:(tf/dt)
% if temp(k,1) < 12.4
%     betah(k,1)= 0;
% elseif temp(k,1) >= 12.4 && temp(k,1)<= 26.1
%    betah(k,1)=0.001044*temp(k,1)*(temp(k,1)-12.286)*sqrt(32.461-temp(k,1));
% else
%     betah(k,1) = 1;
% end
%  end
% betah=0.001044*temp.*(temp-12.286).*sqrt(32.461-temp);%Taghikhani2018
betam=(0.033*temp-0.41);%Lourdes 2015
% betam=-0.9037+0.0729*temp;%Taghikhani2018

mum=(0.8962-0.159*temp+(1.116E-2)*temp.^2-(3.408E-4)*temp.^3+(3.809E-6)*temp.^4);%Putra 2017 e Taghikhani2018
mua=((2.13-0.3797*temp+(2.457E-2)*temp.^2-(6.778E-4)*temp.^3+(6.794E-6)*temp.^4)/7);

% for j=1:1:(tf/dt)
%     if prec(j,1)<50
% %mua(j,1)=2.315-0.419*temp(j,1)+0.02375*(temp(j,1))^2-(7.358E-4)*(temp(j,1))^3+((7.503E-6))*(temp(j,1))^4;%Putra 2017
% mua(j,1)=((2.13-0.3797*temp(j,1)+(2.457E-2)*(temp(j,1))^2-(6.778E-4)*(temp(j,1))^3+(6.794E-6)*(temp(j,1))^4)/7);
%     else
% mua(j,1)=1-exp(-prec(j,1)*log(2)/50);%Putra 2017
%     end
% end
% % plot(t,mua)
alpha=((0.131-(5.723E-2)*temp+(1.164E-2)*temp.^2-(1.341E-3)*temp.^3+(8.723E-5)*temp.^4-(3.017E-6)*temp.^5+(5.153E-8)*temp.^6+(3.42E-10)*temp.^7)/7);
% alpha=-1.847+0.8291*temp-0.1457*temp.^2+(1.305E-2)*temp.^3-(6.461E-4)*temp.^4+(1.796E-5)*temp.^5-(2.61E-7)*temp.^6+(1.551E-9)*temp.^7;%Putra 2017
Cmax=x(3); %ou 1E4
 % Cmax=1E4;
 % Cmax=9.820982413224507E3;

C=((prec./63.2)*Cmax);
for i=1:1:((tf/dt)+1)
f=@(t,y)[(tau^(1-gamma))*(muh*(H-y(1))-(b(i,1)*betah(i,1)*y(1)*y(5))/(H));
    (tau^(1-gamma))*((b(i,1)*betah(i,1)*y(1)*y(5))/(H)-(muh+sigma)*y(2));
    (tau^(1-gamma))*(sigma*y(2)-muh*y(3));
    (tau^(1-gamma))*(alpha(i,1)*y(6)-(b(i,1)*betam(i,1)*y(4)*y(2))/(H)-mum(i,1)*y(4));
    (tau^(1-gamma))*((b(i,1)*betam(i,1)*y(4)*y(2))/(H)-mum(i,1)*y(5));
    (tau^(1-gamma))*((k*delta(i,1)*(1-(y(6)/C(i,1)))*(y(4)+y(5)))-(mua(i,1)+alpha(i,1))*y(6))]; %
size(f(t0, [HS_0; HI_0; HR_0; MS_0; MI_0; A_0]));
%i=i+1;
end

[t,y]=fde12(gamma,f,t0,tf,[HS_0;HI_0;HR_0;MS_0;MI_0;A_0],dt);			% fde12 
% opts = odeset('RelTol',1.e-4);
% 
% [t,y]=ode23s(f,[0 365],[HS_0;HI_0;HR_0;MS_0;MI_0;A_0],opts);
HS_estimado=y(1,:);				% Soluções
HI_estimado=y(2,:);	
HR_estimado=y(3,:);	
MS_estimado=y(4,:);	
MI_estimado=y(5,:);	
A_estimado=y(6,:);	
MA = load('DadosBauru2022.txt');  %dados
n = size(MA);  
t_real=MA(:,1);                                         % primeira coluna (tempo - dias)
HI_real=MA(:,3);   
figure (3)
plot(t_real,HI_real, 'bo')
figure (4)
plot(t,HI_estimado)
figure(2)
plot(t,HI_estimado,'g', t_real, HI_real, 'bo');                   
xlabel('tempo(dias)');                                            
ylabel('População'); 
legend('H_I, com \gamma=1')
title('Cenário 3 - Q_0>1 e R_0>1');
hold on;