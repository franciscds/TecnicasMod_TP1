%---------------------------------------%
%Tecnicas de Modelagem de Sistema      -% 
%       Trabalho Individual 1          -%
%           Exercicio 3                -%
%  Aluno: Francis Carlos dos Santos    -%
%  Matricula: 2012022167               -%
%---------------------------------------%

%Abrir sinal
clear all;
close all;
clc;
matricula = 2012022167;
[xn,fs,nr_samples,rep_start,b1,b2,b3,fsb] = generate_noisy_audio('bee.mp3',matricula);
%prepara objeto de audio a partir de vetor
p = audioplayer(xn,fs);
% play(p);
%  pause(10);
% stop(p);

%% Parte 1 - media coerente a partir de 200 repeticoes
%Numero de repeticoes em um trecho
Ns = 200;
%Numero de amostras em um periodo de 8s
N = nr_samples;
%conjunto de simulacoes Ns cada linha representa uma amostra
for i=1:Ns
y(i,:) = xn((i-1)*N+1:i*N,:);
end
%media de simulacoes ponto a ponto
for i = 1:N
    y_medio(i) = sum(y(:,i))/Ns;
end
% Escuta o Bethoven
% sinal resultante apos a aplicacao da media coerente
 p1 = audioplayer(y_medio,fs);
%  play(p1); 

% uma realizacao do teste
figure
subplot(2,1,1)
plot(1:N,y(1,:),'b-','linewidth',0.5);
xlim([0 N]);
box off;
ylabel('y(t)+e(t)')

% media resultante de Ns testes semelhantes ao mostrado acima
subplot(2,1,2)
plot(1:N,y_medio,'r--','linewidth',0.5);
xlim([0 N]);
box off;
xlabel('t')
ylabel('media[y_i(t)+e_i(t)]')


%% Parte 2 - Media utilizando 20 amostras
for i = 1:N
    y_medio20(i) = sum(y(1:20,i))/20;
end

p2 = audioplayer(y_medio20,fs);
% play(p2);
% pause(10);

% uma realizacao do teste
figure
subplot(2,1,1)
plot(1:N,y(1,:),'b-','linewidth',0.5);
xlim([0 N]);
box off;
ylabel('y(t)+e(t)')

% media resultante de Ns teesta vez stes semelhantes ao mostrado acima
subplot(2,1,2)
plot(1:N,y_medio20,'r--','linewidth',0.5);
xlim([0 N]);
box off;
xlabel('t')
ylabel('Media 20 amostras')
%% Parte 3
% modifica sinal e retira amostras do inicio e final até completar 200
y_mod = zeros(Ns,nr_samples);
for i=1:Ns
	y_mod(i,:) = [y(i,(i+1):N),zeros(1,(N-(N-i)))];
end
% Media Coerente
for i = 1:N
    y_medio_mod(i) = sum(y_mod(:,i))/Ns;
end
% Sinal resultante apos a aplicacao da media coerente
 p3 = audioplayer(y_medio_mod,fs);
% play(p3);
figure
subplot(2,1,1)
plot(1:N,y(1,:),'b-','linewidth',0.5);
xlim([0 N]);
box off;
ylabel('y(t)+e(t)')

% media resultante de Ns testes semelhantes ao mostrado acima
subplot(2,1,2)
plot(1:N,y_medio_mod,'r--','linewidth',0.5);
xlim([0 N]);
box off;
xlabel('t')
ylabel('Media Modificada')