%---------------------------------------%
%Tecnicas de Modelagem de Sistema      -% 
%       Trabalho Individual 1          -%
%           Exercicio 1                -%
%  Aluno: Francis Carlos dos Santos    -%
%  Matricula: 2012022167               -%
%---------------------------------------%
clear all;
% Exercicio 1: Resposta de sistemas e métodos determinísticos
K=1;
tau=1;
h2=1;
h1=1;
K0=1;
sim('system_response_2016_2_parte1_mod')  
R = ScopeData_1.signals(1).values;
Y = ScopeData_1.signals(2).values;
t = ScopeData_1.time';
% degrau de entrada
figure,
subplot(2,1,1)
plot(t,R,'blue','linewidth',2)
box off
xlabel('t')
ylabel('r(t)')
title('Degrau de entrada')
% resposta do sistema ao degrau
subplot(2,1,2)
plot(t,Y,'red','linewidth',2)
box off
title('Saida do sistema')
xlabel('t')
ylabel('y(t)')
%% 1.1 Sistema de primeira Ordem
t0 = 1;
td = 30;
tinf = tout(end);
%amplitude
A = R(end);
%ganho do sistema
K = (Y(end) - Y(t0))./A;
%Calculo de Tau  
Ytau = 0.632*(Y(end)-Y(t0)) + Y(t0);
tauA = find(Y>=Ytau);
tau = tauA(1)/100 -td;

%% 1.2 primeira ordem com atraso puro
tauD = 0.5; % atraso de transorte estimado empiricamente
K2= (Y(end)-Y(t==tauD))/A; % ganho em regime permanente
yd = 0.632*(Y(end)-Y(t==tauD))+Y(t==tauD);
[~,idx2] = min(abs(yd-Y));
tauD = t(idx2); 
tauD = tauD-30;

%% 1.3 Coeficientes sobreamorteciso 
k2 =1;
tm = find(Y >=1);
tm = tm(1)/100 - 30;
% M1
m1 = tm*1/2;

%mi
mi = 1/(tm);
%lambda

lambda = (tm - m1)*mi;
zeta = lambda2zeta(lambda);
eta = lambda2eta(lambda);

wn = (sec(zeta)/sqrt(1 - zeta.^2))*(1/(tm - m1));

h1 = 2*zeta*wn;
h2 = (wn^2);
K0 = K2*(wn^2);
%%

%trazendo degrau para inicio
i = find(t>=30);
ts = t(i)-30; 
u = R(i);
y1 = ScopeData_1.signals(3).values;
y1d= ScopeData_1.signals(4).values;
y2 = ScopeData2.signals.values;
%% Plota grafico comparativo modelo 1 ordem e Y
figure,
plot(t,Y,'m',t,y1,'r','linewidth',2)
box off
xlabel('t')
ylabel('y1(t)')
title('Comparação sistema primeira ordem e saída do sistema ')
%% Plota grafico comparativo modelo 1 ordem com atraso e Y
figure,
plot(t,Y,'m',t,y1d,'b','linewidth',2)
box off
xlabel('t')
ylabel('i1d(t)')
title('Comparação sistema 1 ordem com atraso e saida do sistema ')


%% Plota grafico comparativo modelo 2 ordem sobreamortecido e Y
figure,
plot(t,Y,'blue',t,y2,'red','linewidth',2)
box off
title('Comparação sistema sobreamortecido e saida do sistema')
xlabel('t')
ylabel('y2(t)')
axis tight
%% 1.4 Diagramas no simulink

%% 1.5 Comparar Erro usando mean square error
MSE1_1 = mean((Y-y1).^2);
MSE1_2 = mean((Y-y1d).^2);
MSE1_3 = mean((Y-y2).^2);
%% 1.6
h22=1;
h12=1;
K02=1
sim('system_response_2016_2_parte2') 
y22= ScopeData_2.signals(1).values;
u2 =ScopeData_2.signals(2).values;
t2 = ScopeData_2.time';

figure,
plot(t2,y22,'blue',t2,u2,'red','linewidth',2)
box off
title('Sinal sistema 2 ordem 1.6')
xlabel('t')
ylabel('y(t)')
axis tight

% estimativa de parametros
k22 =1;
tm2 = find(y22 >=1);
tm2 = tm2(1)/100 + 16;
% M1
m12 = tm2*1/2;

%mi
mi2 = 1/(tm2);
%lambda

lambda2 = (tm2 - m12)*mi2;
zeta2 = lambda2zeta(lambda2);
eta2 = lambda2eta(lambda2);

wn2 = (sec(zeta2)/sqrt(1 - zeta2.^2))*(1/(tm2 - m12));

h12 = 2*zeta2*wn2;
h22 = (wn2^2);
K02 = k22*(wn2^2);

% Plota comparação
sim('system_response_2016_2_parte2');
y23 =ScopeData_2.signals(3).values;
figure,
plot(t2,u2,'blue',t2,y23,'red','linewidth',2)
box off
xlabel('t')
ylabel('y(t)')
title('Comparação de Saida gerada e Saida calculada')
legend('y(t) original','y(t) estimado')
%% 1.8 Respondida ---

