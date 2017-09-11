%---------------------------------------%
%Tecnicas de Modelagem de Sistema      -% 
%       Trabalho Individual 1          -%
%           Exercicio 2                -%
%  Aluno: Francis Carlos dos Santos    -%
%  Matricula: 2012022167               -%
%---------------------------------------%
clear all;
close all;
clc;
%% Exercicio 2 - Item A
n = -500:500; % amostras
h(1:500,1) = 0; % h[n] = 0, n < 0 (sistema causal)
u = [zeros(500,1); 1; zeros(500,1)]; % impulso unitario

% resposta h[n] do sistema a um impulso unitario
for i=501:1001
   h(i,1) = (-7/12)*h(i-1,1) - (7/12)*h(i-2,1) + (5/12)*u(i-1,1) + (1/4)*u(i-2,1);
end

figure
stem(n,h); 
box off; 
xlabel('Amostra n'); 
ylabel('h[n]');
legend('h[n]');
title('Resposta ao Impulso do Sistema')

%% Exercicio 2 - Item B

x = randn(501,1); %sendo K = 500
% condicoes iniciais (y[0] = 0, y[1] = 5/12*x[0])
y(1,1) = 0; 
y(2,1) = 5/12*x(1,1); 

% resposta y[n] do sistema ao sinal randomico de entrada
for i = 3:501
    y(i,1) =(-7/12)*y(i-1,1)-(7/12)*y(i-2,1)+(5/12)*x(i-1)+(1/4)*x(i-2);
end

%% Exercicio 2 - Item C

for i = 1:100
    U(i,:) = x(i:(99+i));
end
U = fliplr(U);
% estimativa das 100 primeiras amostras da resposta ao impulso
h_estimado = U\y(100:199);

% compara o valor estimado com o obtido no item 1
figure
stem(1:100,h(501:600),'Color',[0 0.3906 0],'linewidth',2); hold on; 
stem(1:100,h_estimado,'red'); box off; 
xlabel('amostras');
set(gcf,'Position',[250 150 300 350])
legend('h[n] original','h[n] estimado')
title('h[n] original vs h[n] estimado')

%% Exercicio 2 - Item D
y_ruido = y.*(1+0.04.*randn(size(y)));

% estimativa da resposta ao impulso considerando o ruido
h_estimado_ruido = U\y_ruido(100:199);

% compara a estimativa obtida no item anterior com a obtida acima
figure
stem(1:100,h_estimado,'red','linewidth',2); 
hold on 
stem(1:100,h_estimado_ruido,'blue','linewidth',2); 
box off; 
xlabel('amostras');
legend('h[n] Estimado sem ruido','h[n] Estimado com ruido')
title('h[n] Estimado sem ruido VS h[n] Estimado com ruido')

%% Exercicio 3 - Item E
num = [5 3];
den = [12 8 7];
% resposta em frequencia do sistema (H(z) = (5z+3)/(12z^2+7z+7))
figure
freqz(num,den)
title('Resposta em frequencia do sistema')
[h,w] = freqz(num,den);

%% Exercicio 3 - Item F

% ruido branco
r = randn(1024,1);

% condicoes iniciais (yr[0] = 0, yr[1] = 5/12*r[0])
yr(1,1) = 0; 
yr(2,1) = 5/12*r(1,1); 

% resposta y[n] do sistema ao sinal randomico de entrada
for i = 3:1024
    yr(i,1) = (-7/12)*yr(i-1,1)-(7/12)*yr(i-2,1)+(5/12)*r(i-1)+(1/4)*r(i-2);
end
% Caucula trasformada de fourie do sinal com ruido
N=length(yr);
for k=1:N
    FTyr(k)=0;
    for n=1:N
        FTyr(k)=FTyr(k)+yr(n).*exp(-1j.*2.*pi.*(n-1).*(k-1)./N);
    end
end
% Calcula TF de ruido branco
N=length(r);
for k=1:N
    FTr(k)=0;
    for n=1:N
        FTr(k)=FTr(k)+r(n).*exp(-1j.*2.*pi.*(n-1).*(k-1)./N);
    end
end
%Função de transferência
H = FTyr./FTr;
modH = db(H(1:1024/2));
fasH = angle(H(1:1024/2));

for i = 1:length(fasH)
    if (fasH(i) > 2)
        fasH(i) = fasH(i) - 2*pi;
    end
end
fasH = fasH';
figure
subplot(2,1,1)
plot(w,db(h),'blue',w,db(H(1:1024/2)),'red','linewidth',2); box off;
title('Resposta em frequencia do sistema vs resposta em frequencia estimada')
ylabel('Magnitude (dB)')
legend('Magnitude obtida com freqz','Magnitude estimada atraves do ruido branco');

subplot(2,1,2)
plot(w,angle(h),'blue',w,fasH+w,'red','linewidth',2); 
box off;
ylabel('Fase (radianos)')
legend('Fase obtida com freqz', 'Fase estimada atraves do ruido branco');

