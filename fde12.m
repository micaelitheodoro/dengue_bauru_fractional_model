function [t, y] = fde12(alpha,fdefun,t0,tfinal,y0,h)

if nargin < 9
    mu_tol = 1.0e-6 ; 
    if nargin < 8
        mu = 1 ;
        if nargin < 7
            param = [] ;
        end
    end
end

% Verificação da ordem da derivada 
if alpha <= 0
    error('MATLAB:fde12:OrdemNegativa',...
        ['A ordem da derivada tem que ser positiva ' ...
         'ALPHA = %f nao pode ser aceita.'], alpha);
end 

% Verificar o passo do método
if h <= 0
    error('MATLAB:fde12:OrdemNegativa',...
        ['O tamanho do passo tem que ser positivo ' ...
         'H = %e nao pode ser aceito.'], h);
end 

% Armazenar condicoes iniciais
ic.t0 = t0 ;
ic.y0 = y0 ;
ic.m_alpha = ceil(alpha) ;
ic.m_alpha_factorial = factorial(0:ic.m_alpha-1) ;

% Armazenar informacoes do problema
Probl.ic = ic ;
Probl.fdefun = fdefun ;
Probl.problem_size = size(y0,1) ;
Probl.param = param ;

% Verificando o numero de condicoes iniciais do vetor y0
if size(y0,2) < ic.m_alpha
    error('MATLAB:fde12:SemEntradasSuficientes', ...
        ['Nao ha um numero suficiente de entradas. ' ...
        'Ordem ALPHA = %f precisa de %d condicoes iniciais.'], ...
        alpha,ic.m_alpha);
end

% Verificando a compatibilidade do tamanho do problema com o tamanho do
% campo vetorial

f_temp = f_vectorfield(t0,y0(:,1),Probl) ;
if Probl.problem_size ~= size(f_temp,1)
    error('MATLAB:fde12:TamanhoNaoCompativel', ...
        ['O tamanho %d do problema é obtido de acordo com as condicoes iniciais ' ...
         '(i.e. o numero de linhas/colunas de Y0) nao e compativel com ' ...
         'o tamanho %d da saida do vetor FDEFUN. ' ], Probl.problem_size,size(f_temp,1));
end

% Numero dos pontos para avaliar tamanhos e solucoes
r = 16 ;
N = ceil((tfinal-t0)/h) ;
Nr = ceil((N+1)/r)*r ;
Q = ceil(log2(Nr/r)) - 1 ;
NNr = 2^(Q+1)*r ;

% prealocacao
y = zeros(Probl.problem_size,N+1) ;
fy = zeros(Probl.problem_size,N+1) ;
zn_pred = zeros(Probl.problem_size,NNr+1) ;
if mu > 0
    zn_corr = zeros(Probl.problem_size,NNr+1) ;
else
    zn_corr = 0 ;
end

% Avaliando os coeficiente do FracPECE
nvett = 0 : NNr+1 ;
nalpha = nvett.^alpha ;
nalpha1 = nalpha.*nvett ;
PC.bn = nalpha(2:end) - nalpha(1:end-1) ;
PC.an = [ 1 , (nalpha1(1:end-2) - 2*nalpha1(2:end-1) + nalpha1(3:end)) ] ;
PC.a0 = [ 0 , nalpha1(1:end-2)-nalpha(2:end-1).*(nvett(2:end-1)-alpha-1)] ;
PC.halpha1 = h^alpha/gamma(alpha+1) ;
PC.halpha2 = h^alpha/gamma(alpha+2) ;
PC.mu = mu ; PC.mu_tol = mu_tol ;

% inicializando as condicoes iniciais
t = t0 + (0 : N)*h ;
y(:,1) = y0(:,1) ;
fy(:,1) = f_temp ;
[y, fy] = Triangolo(1, r-1, t, y, fy, zn_pred, zn_corr, N, PC, Probl ) ;

% principal processo de iteracao
ff = [0 2 ] ; nx0 = 0 ; ny0 = 0 ;
for q = 0 : Q
    L = 2^q ; 
    [y, fy] = DisegnaBlocchi(L, ff, r, Nr, nx0+L*r, ny0, t, y, fy, ...
                             zn_pred, zn_corr, N, PC, Probl ) ;
    ff = [ff ff]  ; ff(end) = 4*L ;
end

% Avaliacao da solucao no t_final
if tfinal < t(N+1) 
    c = (tfinal - t(N))/h ;
    t(N+1) = tfinal ;
    y(:,N+1) = (1-c)*y(:,N) + c*y(:,N+1) ;
end
t = t(1:N+1) ; y = y(:,1:N+1) ;

end

% r : dimensão of do triangulo ou quadrado basico
% L : fator de redimensionalizacao dos quadrados
function [y, fy] = DisegnaBlocchi(L, ff, r, Nr, nx0, ny0, t, y, fy, ...
                                  zn_pred, zn_corr, N , PC, Probl)

nxi = nx0 ; nxf = nx0 + L*r - 1 ;
nyi = ny0 ; nyf = ny0 + L*r - 1 ;
is = 1 ;
s_nxi(is) = nxi ; s_nxf(is) = nxf ; s_nyi(is) = nyi ; s_nyf(is) = nyf ;

i_triangolo = 0 ; stop = 0 ;
while ~stop
    
    stop = nxi+r-1 == nx0+L*r-1 | (nxi+r-1>=Nr-1) ; % para quando o triângulo atual termina no final do quadrado
    
    [zn_pred, zn_corr] = Quadrato(nxi, nxf, nyi, nyf, fy, zn_pred, zn_corr, PC, Probl) ;
    
    [y, fy] = Triangolo(nxi, nxi+r-1, t, y, fy, zn_pred, zn_corr, N, PC, Probl) ;
    i_triangolo = i_triangolo + 1 ;
    
    if ~stop
        if nxi+r-1 == nxf   %O triângulo termina onde termina o quadrado -> nivela para baixo
            i_Delta = ff(i_triangolo) ;
            Delta = i_Delta*r ;
            nxi = s_nxf(is)+1 ; nxf = s_nxf(is)  + Delta ;
            nyi = s_nxf(is) - Delta +1; nyf = s_nxf(is)  ;
            s_nxi(is) = nxi ; s_nxf(is) = nxf ; s_nyi(is) = nyi ; s_nyf(is) = nyf ;
        else % O triângulo termina antes do quadrado -> faz um quadrado próximo a ele
            nxi = nxi + r ; nxf = nxi + r - 1 ; nyi = nyf + 1 ; nyf = nyf + r  ;
            is = is + 1 ;
            s_nxi(is) = nxi ; s_nxf(is) = nxf ; s_nyi(is) = nyi ; s_nyf(is) = nyf ;
        end
    end
    
end

end

function [zn_pred, zn_corr] = Quadrato(nxi, nxf, nyi, nyf, fy, zn_pred, zn_corr, PC, Probl)

coef_beg = nxi-nyf ; coef_end = nxf-nyi+1 ;
funz_beg = nyi+1 ; funz_end = nyf+1 ;

%Avaliacao do segmento da convolucao
vett_coef = PC.bn(coef_beg:coef_end) ;
vett_funz = [fy(:,funz_beg:funz_end) , zeros(Probl.problem_size,funz_end-funz_beg+1) ] ;
zzn_pred = real(FastConv(vett_coef,vett_funz)) ;
zn_pred(:,nxi+1:nxf+1) = zn_pred(:,nxi+1:nxf+1) + zzn_pred(:,nxf-nyf+1-1:end-1) ; 


if PC.mu > 0
    vett_coef = PC.an(coef_beg:coef_end) ;
    if nyi == 0 % avaliando o menor quadrado
        vett_funz = [zeros(Probl.problem_size,1) , fy(:,funz_beg+1:funz_end) , zeros(Probl.problem_size,funz_end-funz_beg+1) ] ;
    else % avaliando todos menos o menor
        vett_funz = [ fy(:,funz_beg:funz_end) , zeros(Probl.problem_size,funz_end-funz_beg+1) ] ;
    end
    zzn_corr = real(FastConv(vett_coef,vett_funz)) ;
    zn_corr(:,nxi+1:nxf+1) = zn_corr(:,nxi+1:nxf+1) + zzn_corr(:,nxf-nyf+1:end) ;
else
    zn_corr = 0 ;
end

end

function [y, fy] = Triangolo(nxi, nxf, t, y, fy, zn_pred, zn_corr, N, PC, Probl)

for n = nxi : min(N,nxf)
    
    % avaliando o preditor
    Phi = zeros(Probl.problem_size,1) ;
    if nxi == 1 % Caso do primeiro triangulo
        j_beg = 0 ;
    else % caso dos demais triangulos sem ser o primeiro
        j_beg = nxi ;
    end
    for j = j_beg : n-1
        Phi = Phi + PC.bn(n-j)*fy(:,j+1) ;
    end
    
    St = StartingTerm(t(n+1),Probl.ic) ;
    y_pred = St + PC.halpha1*(zn_pred(:,n+1) + Phi) ;
    f_pred = f_vectorfield(t(n+1),y_pred,Probl) ; 
    
    if PC.mu == 0
        y(:,n+1) = y_pred ;
        fy(:,n+1) = f_pred ;
    else
        j_beg = nxi ;
        Phi = zeros(Probl.problem_size,1) ;
        for j = j_beg : n-1
            Phi = Phi + PC.an(n-j+1)*fy(:,j+1) ;
        end
        Phi_n = St + PC.halpha2*(PC.a0(n+1)*fy(:,1) + zn_corr(:,n+1) + Phi) ;
        yn0 = y_pred ; fn0 = f_pred ;
        stop = 0 ; mu_it = 0 ;
        while ~stop
            yn1 = Phi_n + PC.halpha2*fn0 ;
            mu_it = mu_it + 1 ;
            if PC.mu == Inf
                stop = norm(yn1-yn0,inf) < PC.mu_tol ;
                if mu_it > 100 && ~stop
                    warning('MATLAB:fde12:NaoConvergencia',...
                        strcat('Foi solicitada a execução de iterações de corretor até a convergência, mas ', ...
                        'o processo nao converge para a tolerancia %e em 100 iteracoes'),PC.mu_tol) ;
                    stop = 1 ;
                end
            else
                stop = mu_it == PC.mu ;
            end
            fn1 = f_vectorfield(t(n+1),yn1,Probl) ;             
            yn0 = yn1 ; fn0 = fn1 ;
        end
        y(:,n+1) = yn1 ;
        fy(:,n+1) = fn1 ;
    end
end

end

function z = FastConv(x, y)

r = length(x) ; problem_size = size(y,1) ;

z = zeros(problem_size,r) ; 
X = fft(x,r) ;
for i = 1 : problem_size
    Y = fft(y(i,:),r) ;
    Z = X.*Y ;
    z(i,:) = ifft(Z,r) ;
end

end

function f = f_vectorfield(t,y,Probl)

if isempty(Probl.param)
    f = feval(Probl.fdefun,t,y) ;
else
    f = feval(Probl.fdefun,t,y,Probl.param) ;
end

end

function ys = StartingTerm(t,ic)

ys = zeros(size(ic.y0,1),1) ;
for k = 1 : ic.m_alpha
    ys = ys + (t-ic.t0)^(k-1)/ic.m_alpha_factorial(k)*ic.y0(:,k) ;
end

end