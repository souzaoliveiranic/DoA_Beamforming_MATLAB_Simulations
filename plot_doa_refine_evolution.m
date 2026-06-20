% =========================================================================
% GRAFICO: evolucao do erro de DoA por transmissao durante o refino C<->DoA.
% -------------------------------------------------------------------------
% Dentro de UMA execucao multishot, mostra como o erro de direcao de cada
% transmissao evolui ao longo das iteracoes de refino. Uma curva por
% transmissao (1a, 2a, ...). Eixo x = iteracao de refino (0 = inicial).
%
% Para a evolucao ser visivel, parte-se de direcoes iniciais COM ERRO
% (perturbadas), simulando um estimador inicial pobre; o refino conjunto
% C<->DoA deve reduzir o erro de cada transmissao ao longo das iteracoes.
%
% AUTONOMO. Depende de: utils.steering_vec_uca, utils.fsk2_mod, doa_kw_uca,
%   estimate_C_circulant_pool.
% Saida: doa_refine_evolution.png
% =========================================================================
clear; close all; clc;

M=8; fc=500e6; c0=3e8; lambda=c0/fc; radius=0.2*lambda;
K_kw=2100; Rs=9600; sps=30; alpha=0.3; span=8; fd=4.8e3;
SNR_dB=10;                 % SNR alto p/ isolar a dinamica do refino (pouco ruido)
sc_lambda=1e-2;
L=6;                       % numero de transmissoes (curvas)
NREF=8;                    % iteracoes de refino a mostrar
init_err_deg=10;           % erro de direcao inicial (perturbacao)
use_safeguard=true;        % true: refino com criterio de parada
seed=20260622;             % reprodutivel

cc=[1,0.40*exp(1j*pi/6),0.15,0.05,0.02,0.05,0.15,0.40*exp(-1j*pi/6)].';
cc(M/2+1)=real(cc(M/2+1));
C_true=zeros(M); for col=1:M,C_true(:,col)=circshift(cc,col-1);end

% preambulo fixo
rng(20260101); Nsym=ceil(K_kw/sps)+span+4;
qv=utils.fsk2_mod(Nsym,Rs,sps,alpha,span,fd); qv=qv(:); qv=qv(1:K_kw);
qv=qv/sqrt(mean(abs(qv).^2));
beta=2*pi*(0:M-1).'/M;

% --- gera L transmissoes (sem interferidor) ---
rng(seed);
phis_true=zeros(1,L); Xw=cell(1,L); B=zeros(M,L);
for l=1:L
  phi_sig=randi([-180,179]); phis_true(l)=phi_sig;
  a_s=utils.steering_vec_uca(M,radius,lambda,90,phi_sig);
  Xn=sqrt(10^(-SNR_dB/10)/2)*(randn(M,K_kw)+1j*randn(M,K_kw));
  Xw{l}=C_true*(a_s*qv.')+Xn;
  B(:,l)=Xw{l}*conj(qv)/(qv'*qv);
end

% --- direcao inicial COM ERRO (it 0) ---
phi=phis_true+init_err_deg*(2*rand(1,L)-1);

% historico de erro por transmissao x iteracao
err_hist=zeros(L,NREF+1);
A=zeros(M,L); for l=1:L,A(:,l)=utils.steering_vec_uca(M,radius,lambda,90,phi(l));end
C=estimate_C_circulant_pool(B,A,M,sc_lambda,eye(M));
ee=abs(phi-phis_true);ee=min(ee,360-ee); err_hist(:,1)=ee(:);
% residuo inicial
rescur=0;nb=0;
for l=1:L,Ca=C*A(:,l);be=(Ca'*B(:,l))/(Ca'*Ca+eps);e=B(:,l)-be*Ca;rescur=rescur+real(e'*e);nb=nb+real(B(:,l)'*B(:,l));end
rescur=rescur/(nb+eps);

for it=1:NREF
  D=(C'*C+1e-9*eye(M))\C'; Aref=zeros(M,L); phiref=zeros(1,L);
  for l=1:L,[~,pk,~]=doa_kw_uca(D*Xw{l},qv.',radius,lambda,beta);phiref(l)=pk(1);Aref(:,l)=utils.steering_vec_uca(M,radius,lambda,90,phiref(l));end
  Ctry=estimate_C_circulant_pool(B,Aref,M,sc_lambda,eye(M));
  rtry=0;nb=0;
  for l=1:L,Ca=Ctry*Aref(:,l);be=(Ca'*B(:,l))/(Ca'*Ca+eps);e=B(:,l)-be*Ca;rtry=rtry+real(e'*e);nb=nb+real(B(:,l)'*B(:,l));end
  rtry=rtry/(nb+eps);
  accept = (~use_safeguard) || (rtry<rescur*(1-1e-4));
  if accept, C=Ctry; phi=phiref; rescur=rtry; end
  ee=abs(phi-phis_true);ee=min(ee,360-ee); err_hist(:,it+1)=ee(:);
end

% --- plot: uma curva por transmissao ---
f=figure('Name','Evolucao do erro de DoA no refino','NumberTitle','off','Position',[60 60 760 500]);
hold on; grid on;
cmap=lines(L);
leg=cell(1,L);
for l=1:L
  plot(0:NREF, err_hist(l,:), '-o', 'Color',cmap(l,:), 'LineWidth',1.6, 'MarkerSize',5,'MarkerFaceColor',cmap(l,:));
  leg{l}=sprintf('transmissao %d (\\phi=%d^o)', l, phis_true(l));
end
xlabel('Iteracao de refino'); ylabel('Erro de DoA |\phi^{hat}-\phi^{true}| (graus)');
title(sprintf('Evolucao do erro de DoA por transmissao durante o refino  |  L=%d, SNR=%+d', L, SNR_dB));
legend(leg, 'Location','northeast','FontSize',8);
xlim([-0.3 NREF+0.3]); ylim([0 max(err_hist(:))*1.1+0.5]);
exportgraphics(f, fullfile(pwd,'doa_refine_evolution.png'), 'Resolution', 200);
fprintf('-> doa_refine_evolution.png\n');
fprintf('Erro medio inicial=%.2f  final=%.2f deg\n', mean(err_hist(:,1)), mean(err_hist(:,end)));
