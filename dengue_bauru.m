%Modelo Dengue - Fracionário
clc
clear all
dt=0.001;		
tf =365;	    
t=0:dt:tf;		
t0=0;

temp=((11.19/2)*cos((2*pi*t)/365)+22.77)'; %Função de temperatura ao longo do ano                      
prec=((63.2/2)*cos((2*pi*t)/365)+(63.2/2))';%Função de precipitação ao longo do ano
HS_0=379138;
HI_0=8;
HR_0=x(2);
MS_0=568719;
MI_0=0;
A_0=0;
muh=3.9E-5;
sigma=(1/7);
%H=379146; %População Bauru
H = ones(3650001,1)*379146;
k=0.8;
delta=-15.837+(1.2897*temp)-0.0163*(temp).^2;
plot(t,prec)
%temp.^2
b=(0.056*delta);   

gamma=1;
tau=1;
betah=0.4;
betam=(0.033*temp-0.41);

mum=(0.8962-0.159*temp+(1.116E-2)*temp.^2-(3.408E-4)*temp.^3+(3.809E-6)*temp.^4);

mua=((2.13-0.3797*temp+(2.457E-2)*temp.^2-(6.778E-4)*temp.^3+(6.794E-6)*temp.^4)/7);

alpha=((0.131-(5.723E-2)*temp+(1.164E-2)*temp.^2-(1.341E-3)*temp.^3+(8.723E-5)*temp.^4-(3.017E-6)*temp.^5+(5.153E-8)*temp.^6+(3.42E-10)*temp.^7)/7);

Cmax=1E6;
C=((prec/63.2)*Cmax);

f=@(t,y)[(tau^(1-gamma))*(muh*(H(fix(round(t,3)/dt+1))-y(1))-((b(fix(round(t,3)/dt+1))*betah).*y(1).*y(5))./(H(fix(round(t,3)/dt+1))));
   (tau^(1-gamma))*(((b(fix(round(t,3)/dt+1))*betah).*y(1).*y(5))./(H(fix(round(t,3)/dt+1)))-(muh+sigma)*y(2));
    (tau^(1-gamma))*(sigma*y(2)-muh*y(3));
    (tau^(1-gamma))*(alpha(fix(round(t,3)/dt+1)).*y(6)-(b(fix(round(t,3)/dt+1)).*betam(fix(round(t,3)/dt+1)).*y(4).*y(2))./(H(fix(round(t,3)/dt+1)))-mum(fix(round(t,3)/dt+1)).*y(4));
    (tau^(1-gamma))*(((b(fix(round(t,3)/dt+1)).*betam(fix(round(t,3)/dt+1)).*y(4).*y(2))./(H(fix(round(t,3)/dt+1))))-mum(fix(round(t,3)/dt+1)).*y(5));
    (tau^(1-gamma))*((((k*delta(fix(round(t,3)/dt+1))).*(1-(y(6)./C(fix(round(t,3)/dt+1))))).*(y(4)+y(5)))-(mua(fix(round(t,3)/dt+1))+alpha(fix(round(t,3)/dt+1))).*y(6))]; %
size(f(t0, [HS_0; HI_0; HR_0; MS_0; MI_0; A_0]))

[f,y]=fde12(gamma,f,t0,tf,[HS_0;HI_0;HR_0;MS_0;MI_0;A_0],dt);			% fde12 
% opts = odeset('RelTol',1.e-4);
% 
% [t,y]=ode23s(f,[0 365],[HS_0;HI_0;HR_0;MS_0;MI_0;A_0],opts);
HS_estimado=y(1,:);				% Soluções
HI_estimado=y(2,:);	
HR_estimado=y(3,:);	
MS_estimado=y(4,:);	
MI_estimado=y(5,:);	
A_estimado=y(6,:);	

% MA = load('DadosBauru2022.txt');  %colocar o csv
% n = size(MA);  
% t_real=MA(:,1);                                         % primeira coluna (tempo - dias)
% HI_real=MA(:,3);   
% figure(2)
% plot(t,HI_estimado,'g', t_real, HI_real, 'bo');                   
% xlabel('tempo(dias)');                                            
% ylabel('População'); 
% legend('H_I, com \gamma=1')
% title('Cenário 3 - Q_0>1 e R_0>1');
% hold on;