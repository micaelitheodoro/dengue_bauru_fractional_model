function diff = fit_simp(x,MA)			
t_real=MA(:,1);						% primeira coluna (tempo - dias)
HI_real=MA(:,3);						% terceira coluna (Y - casos)
dengue_bauru
n = size(MA);
for i=1:n(1)
  diff(i) = sqrt(((HI_real(i)-HI_estimado(round(t_real(i)/dt)+1))^2)); % Residuo, retirei a comparação com D_real, mas basta somar com a diferença quadrática
end
end
