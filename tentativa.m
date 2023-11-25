%Modelo Dengue - Fracionário
%clc
%clear all
dt=0.0001;		
tf =315;	    
t=0:dt:tf;		
t0=0;
%temp=((11.19/2)*cos((2*pi*t)/365)+22.77)'; %Função de temperatura ao longo do ano 365
temp=((11.45/2)*cos((2*pi*t)/365)+22.33)'; %Teste com 315 dias contando a partir do dia 50
% a0=22.8488;%    22.6783    23.0194
% a1=2.8910;%     2.6336     3.1483 
% b1=1.8890;%     1.5701     2.2078 
% a2=-0.4864;%    -0.7173    -0.2556
% b2=0.6046;%     0.3693     0.8399 
% a3=-0.6866;%    -0.9762    -0.3970
% b3=-0.7965;%    -1.0716    -0.5215
% a4=-0.0889;%    -0.4100    0.2321 
% b4=0.5964;%     0.3591     0.8338 
% a5=-0.0210 ;%   -0.2548    0.2127 
% b5=-0.2285;%    -0.4631    0.0060 
% a6=0.2138;%     -0.0194    0.4470 
% b6=0.0481;%     -0.1936    0.2898 
% a7=0.2763 ;%    0.0020     0.5506 
% b7=0.0003 ;%    -0.2504    0.2510 
% a8=-0.7948;   % -1.0609    -0.5287
% b8=0.2863 ;%    -0.2285    0.8010 
% w=0.0177 ;%    0.0173     0.0181 
% temp=(a0+a1*cos(t*w)+b1*sin(t*w)+a2*cos(2*t*w)+b2*sin(2*t*w)+a3*cos(3*t*w)+b3*sin(3*t*w)+a4*cos(4*t*w)+b4*sin(4*t*w)+a5*cos(5*t*w)+b5*sin(5*t*w)+a6*cos(6*t*w)+b6*sin(6*t*w)+a7*cos(7*t*w)+b7*sin(7*t*w)+a8*cos(8*t*w)+b8*sin(8*t*w))';
% plot(t,temp);
% c0=1.0267e11 ;%    -4.0826e+15    4.0828e+15
% c1=-1.3689e11  ;%  -5.4438e+15    5.4435e+15
% d1=5.1849e08     ;%-1.5463e+13    1.5464e+13
% c2=3.4222e10     ;%-1.3609e+15    1.3609e+15
% d2=-2.5924e+8 ;%   -7.7320e+12    7.7315e+12
% w1=-2.0180e-5 ;%   -0.2006        0.2006    

%prec=(c0 + c1*cos(t*w1) + d1*sin(t*w1) + c2*cos(2*t*w1) + d2*sin(2*t*w1))';
%plot(t,prec)
prec=((63.2/2)*cos((2*pi*t)/365)+(63.2/2))';%Função de precipitação ao longo do ano
HS_0=379110;
%HS_0=379138;%para 365 dias
HI_0=36;
%HI_0=8; %para 365 dias
HR_0=0;%HR_0=0; %para 365 dias
MS_0=568588;
MI_0=54;
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
%betah=0.4;
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
% alpha=((0.131-(5.723E-2)*temp+(1.164E-2)*temp.^2-(1.341E-3)*temp.^3+(8.723E-5)*temp.^4-(3.017E-6)*temp.^5+(5.153E-8)*temp.^6+(3.42E-10)*temp.^7)/7);
alpha=-1.847+0.8291*temp-0.1457*temp.^2+(1.305E-2)*temp.^3-(6.461E-4)*temp.^4+(1.796E-5)*temp.^5-(2.61E-7)*temp.^6+(1.551E-9)*temp.^7;%Putra 2017
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
MA = load('DadosBauru2022_315.csv');  %dados
n = size(MA);  
t_real=MA(:,1);                                         % primeira coluna (tempo - dias)
HI_real=MA(:,3);   


MATRIZ_ICC_CASOS=zeros(53,2);
%MATRIZ_ICC_MORTES=zeros(61,2);
for i=1:n(1)
MATRIZ_ICC_CASOS(i,1)=HI_real(i);
MATRIZ_ICC_CASOS(i,2)=HI_estimado(round(t_real(i)/dt)+1);
%MATRIZ_ICC_MORTES(i,1)=D_real(i);
%MATRIZ_ICC_MORTES(i,2)=D_estimado(round(t_real(i)/dt)+1);
end
ICC_CASOS=ICC(MATRIZ_ICC_CASOS,'C-1')
%ICC_MORTES=ICC(MATRIZ_ICC_MORTES,'C-1')
% 
% figure (3)
% plot(t_real,HI_real, 'bo')
% figure (4)
% plot(t,HI_estimado)
% figure(2)
% plot(t,HI_estimado,'g', t_real, HI_real, 'bo');                   
% xlabel('tempo(dias)');                                            
% ylabel('População'); 
% legend('H_I, com \gamma=1')
% title('Cenário 3 - Q_0>1 e R_0>1');
% hold on;