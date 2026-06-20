% =========================================================================
% TESTE: importancia do criterio de parada (salvaguarda) no refino conjunto.
% -------------------------------------------------------------------------
% Compara refino COM salvaguarda vs SEM, medindo Frobenius e RMSE de DoA a
% cada iteracao de refino, em DOIS cenarios:
%   (A) direcao inicial BOA  (KW cru ~1 grau)  -> refino tem pouco a corrigir
%   (B) direcao inicial RUIM (erro ~8 graus)   -> refino precisa corrigir
%
% Responde: SEM salvaguarda, o refino MANTEM o valor bom ou PIORA?
%   - Cenario A mostra o RISCO (degradar uma solucao ja boa).
%   - Cenario B mostra o GANHO (consertar uma solucao ruim) e se o refino
%     cego "passa do ponto".
%
% AUTONOMO. Usa o pool diretamente p/ controlar a direcao inicial (Cenario B);
% Cenario A usa a funcao multishot completa.
% Depende de: utils.steering_vec_uca, utils.fsk2_mod, doa_kw_uca,
%   estimate_C_circulant_pool.
% =========================================================================
clear; close all; clc;

M=8; fc=500e6; c0=3e8; lambda=c0/fc; radius=0.2*lambda;
K_kw=2100; Rs=9600; sps=30; alpha=0.3; span=8; fd=4.8e3;
SNR_dB=3; ISR_dB=3; sc_lambda=1e-2; min_sep=25; L=8;
NREF=4; n_rep=40; rho=0.0;       % cenario limpo (a salvaguarda e' sobre DoA)
init_err_deg=8;                  % erro de direcao inicial no Cenario B

cc=[1,0.40*exp(1j*pi/6),0.15,0.05,0.02,0.05,0.15,0.40*exp(-1j*pi/6)].';
cc(M/2+1)=real(cc(M/2+1));
C_true=zeros(M); for col=1:M,C_true(:,col)=circshift(cc,col-1);end
CtN=C_true/C_true(1,1); froT=norm(CtN,'fro');

rng(20260101); Nsym=ceil(K_kw/sps)+span+4;
qv=utils.fsk2_mod(Nsym,Rs,sps,alpha,span,fd); qv=qv(:); qv=qv(1:K_kw);
qv=qv/sqrt(mean(abs(qv).^2)); rng('shuffle');
beta=2*pi*(0:M-1).'/M;

% acumuladores [modo(1=safe,2=nosafe)][nref+1], por cenario
froA=zeros(2,NREF+1); doaA=zeros(2,NREF+1);
froB=zeros(2,NREF+1); doaB=zeros(2,NREF+1);

for rep=1:n_rep
 % gera L transmissoes
 Xw=cell(1,L); phis_true=zeros(1,L); B=zeros(M,L);
 for l=1:L
  phi_sig=randi([-180,179]); phis_true(l)=phi_sig;
  pim=randi([-180,179]); while min(abs(pim-phi_sig),360-abs(pim-phi_sig))<min_sep, pim=randi([-180,179]); end
  a_s=utils.steering_vec_uca(M,radius,lambda,90,phi_sig);
  qn=qv.'/sqrt(mean(abs(qv).^2)); w=(randn(1,K_kw)+1j*randn(1,K_kw))/sqrt(2);
  vint=rho*qn+sqrt(1-rho^2)*w; vint=vint/sqrt(mean(abs(vint).^2));
  a_i=utils.steering_vec_uca(M,radius,lambda,90,pim);
  Xw{l}=C_true*(a_s*qv.'+a_i*(10^(ISR_dB/20)*vint))+sqrt(10^(-SNR_dB/10)/2)*(randn(M,K_kw)+1j*randn(M,K_kw));
  B(:,l)=Xw{l}*conj(qv)/(qv'*qv);
 end

 % ---- direcao inicial dos dois cenarios ----
 phi0_good=zeros(1,L);
 for l=1:L,[~,pk,~]=doa_kw_uca(Xw{l},qv.',radius,lambda,beta);phi0_good(l)=pk(1);end
 phi0_bad=phis_true+init_err_deg*(2*rand(1,L)-1);

 for sc=1:2     % 1=A(bom), 2=B(ruim)
  if sc==1, phi0=phi0_good; else, phi0=phi0_bad; end
  for mode=1:2  % 1=safe, 2=nosafe
   phi=phi0; A=zeros(M,L); for l=1:L,A(:,l)=utils.steering_vec_uca(M,radius,lambda,90,phi(l));end
   C=estimate_C_circulant_pool(B,A,M,sc_lambda,eye(M));
   rescur=0; nb=0;
   for l=1:L, Ca=C*A(:,l); be=(Ca'*B(:,l))/(Ca'*Ca+eps); e=B(:,l)-be*Ca; rescur=rescur+real(e'*e); nb=nb+real(B(:,l)'*B(:,l)); end
   rescur=rescur/(nb+eps);
   % registra nref=0
   ee=abs(phi-phis_true);ee=min(ee,360-ee);
   if sc==1, froA(mode,1)=froA(mode,1)+norm(C/C(1,1)-CtN,'fro')/froT; doaA(mode,1)=doaA(mode,1)+mean(ee.^2);
   else,     froB(mode,1)=froB(mode,1)+norm(C/C(1,1)-CtN,'fro')/froT; doaB(mode,1)=doaB(mode,1)+mean(ee.^2); end
   for it=1:NREF
    D=(C'*C+1e-9*eye(M))\C'; Aref=zeros(M,L); phiref=zeros(1,L);
    for l=1:L,[~,pk,~]=doa_kw_uca(D*Xw{l},qv.',radius,lambda,beta);phiref(l)=pk(1);Aref(:,l)=utils.steering_vec_uca(M,radius,lambda,90,phiref(l));end
    Ctry=estimate_C_circulant_pool(B,Aref,M,sc_lambda,eye(M));
    rtry=0; nb=0;
    for l=1:L, Ca=Ctry*Aref(:,l); be=(Ca'*B(:,l))/(Ca'*Ca+eps); e=B(:,l)-be*Ca; rtry=rtry+real(e'*e); nb=nb+real(B(:,l)'*B(:,l)); end
    rtry=rtry/(nb+eps);
    accept = (mode==2) || (rtry<rescur*(1-1e-4));   % mode2=sem salvaguarda: aceita sempre
    if accept, C=Ctry; phi=phiref; rescur=rtry; end
    ee=abs(phi-phis_true);ee=min(ee,360-ee);
    if sc==1, froA(mode,it+1)=froA(mode,it+1)+norm(C/C(1,1)-CtN,'fro')/froT; doaA(mode,it+1)=doaA(mode,it+1)+mean(ee.^2);
    else,     froB(mode,it+1)=froB(mode,it+1)+norm(C/C(1,1)-CtN,'fro')/froT; doaB(mode,it+1)=doaB(mode,it+1)+mean(ee.^2); end
   end
  end
 end
end
froA=froA/n_rep; doaA=sqrt(doaA/n_rep); froB=froB/n_rep; doaB=sqrt(doaB/n_rep);

pr=@(v) fprintf('%6.3f ',v); prd=@(v) fprintf('%6.2f ',v);
fprintf('\n===== Cenario A: direcao inicial BOA (KW cru) =====\n');
fprintf(' nref:               '); fprintf('%6d ',0:NREF); fprintf('\n');
fprintf(' Frob COM salvaguarda:'); pr(froA(1,:)); fprintf('\n');
fprintf(' Frob SEM salvaguarda:'); pr(froA(2,:)); fprintf('\n');
fprintf(' RMSE COM salvaguarda:'); prd(doaA(1,:)); fprintf('\n');
fprintf(' RMSE SEM salvaguarda:'); prd(doaA(2,:)); fprintf('\n');
fprintf('\n===== Cenario B: direcao inicial RUIM (~%ddeg) =====\n',init_err_deg);
fprintf(' nref:               '); fprintf('%6d ',0:NREF); fprintf('\n');
fprintf(' Frob COM salvaguarda:'); pr(froB(1,:)); fprintf('\n');
fprintf(' Frob SEM salvaguarda:'); pr(froB(2,:)); fprintf('\n');
fprintf(' RMSE COM salvaguarda:'); prd(doaB(1,:)); fprintf('\n');
fprintf(' RMSE SEM salvaguarda:'); prd(doaB(2,:)); fprintf('\n');
