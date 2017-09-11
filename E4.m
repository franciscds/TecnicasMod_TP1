%---------------------------------------%
%Tecnicas de Modelagem de Sistema      -% 
%       Trabalho Individual 1          -%
%           Exercicio 3                -%
%  Aluno: Francis Carlos dos Santos    -%
%  Matricula: 2012022167               -%
%---------------------------------------%
matricula = 2012022167;
[xn,fs,nr_samples,rep_start,b1,b2,b3,fsb] = generate_noisy_audio('bee.mp3',matricula);
%conjunto de simulacoes Ns cada linha representa uma amostra
Ns = 200;
N = nr_samples;
for i=1:Ns
y(i,:) = xn((i-1)*N+1:i*N,:);
end
%media de simulacoes ponto a ponto
for i = 1:N
    y_medio(i) = sum(y(:,i))/Ns;
end
b0 = y_medio;
p1 = audioplayer(b1,fsb);
% play(p1);
% pause(10);
% stop(p1)
p2 = audioplayer(b2,14700);
% play(p2);
% pause(10);
stop(p2)
p3 = audioplayer(b3,7350);
% play(p3);
% pause(10);
% stop(p3)
%% Item 2 - Estimar Atraso
FS = fsb;
%Temos 3 sinais com frequencia de amostragems multiplas
fb = 7350;

b1_mod = b1(1:FS/fb:3*N);
b2_mod = b2(1:14700/fb:2*N);
b3_mod = b3(1:N);
s1 = [b1_mod b0'];
s2 = [b2_mod b0'];
s3 = [b3_mod b0'];
C_b1b0=xcov(s1,'unbiased'); % FCC entre os sinais b1 e b0
C_b2b0=xcov(s2,'unbiased'); % FCC entre os sinais b2 e b0
C_b3b0=xcov(s3,'unbiased'); % FCC entre os sinais b3 e b0

plot(C_b1b0,'b-') 
hold on
plot(C_b2b0,'r-')
plot(C_b3b0,'m-')
xlabel('Amostras')
legend('b1b0','b2b0','b3b0')
