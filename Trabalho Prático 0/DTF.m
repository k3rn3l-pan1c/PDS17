%% Trabalho Prático 0 - Exercícios de revisão DFT
clear; clc; close all;

%% Ex 1
% [ANALITICO]

%% Ex 2
% Número de amostras e vetor de indices
N = 2;
n = 0:N-1;

% Sinal 'a', 'b' e respetivas DFTs
a = [1 0];
b = [1 1];

A = fft(a);
B = fft(b);


figure(1)
subplot(221)
plot(n, a, 'x', 'MarkerSize', 10, 'LineWidth', 4);
title('Sinal a')
xlabel('Indíce_{amostra}')
ylabel('Amplitude')

subplot(222)
plot(n, A, 'x', 'MarkerSize', 10, 'LineWidth', 4)
title('DFT sinal a')
xlabel('Indíce_{amostra}')
ylabel('Amplitude')

subplot(223)
plot(n, b, 'x', 'MarkerSize', 10, 'LineWidth', 4)
title('Sinal b')
xlabel('Indíce_{amostra}')
ylabel('Amplitude')

subplot(224)
plot(n, B, 'x', 'MarkerSize', 10, 'LineWidth', 4)
title('DFT sinal b')
xlabel('Indíce_{amostra}')
ylabel('Amplitude')

%% Exercício 3
% Número de amostras e vetor de indices
N = 4;
n = 0:N-1;

x = [1 1 0 1];
X = fft(x);

figure(2)
subplot(211)
stem(n, x, 'x', 'MarkerSize', 10, 'LineWidth', 3);
title('Sinal a')
xlabel('Indíce_{amostra}')
ylabel('Amplitude')
grid on

subplot(212)
stem(n, X, 'x', 'MarkerSize', 10, 'LineWidth', 3)
title('DFT sinal a')
xlabel('Indíce_{amostra}')
ylabel('Amplitude')
grid on
%% Exercício 4
% Número de amostras e vetor de indices
N = 100;
n = 0:N-1;

x = cos(2*pi*40 * n ./ N);


figure(3)
stem(n, x, '.', 'MarkerSize', 20)
hold on
plot(n, x, '--')
hold off
title('Sinal x');
xlabel('Amostras');
ylabel('Amplitude');
legend('Valor amostras', 'sinal interpolado')

%% Alínea a)
X = fft(x);

figure(4)
subplot(121)
stem(n, abs(X), '.', 'MarkerSize', 20)
title('DFT(x) centrada em Fa/2');
xlabel('Amostras');
ylabel('Amplitude');

subplot(122)
stem(-N/2:N/2-1, abs(fftshift(X)), '.', 'MarkerSize', 20)
title('DFT(x) centrada em 0');
xlabel('Amostras');
ylabel('Amplitude');

%% Alínea b)
% Número de amostras por símbolo e vetor de indices
L = 2;
l = 0:2*N-1;

% Upsampling (L amostras por símbolo)
z = zeros(1, N*L);
z(1:L:end) = x;

Z = fft(z);

figure(5)
subplot(2,2,[1,2])
stem(l(1:L:end), x, 'x')
hold on
stem(l, z, 's')
hold off
title('Sinal z vs x')
xlabel('Indíce das amostras')
ylabel('Amplitude')

subplot(223)
stem(l, abs(Z), 'x', 'MarkerSize', 10)
hold on
stem(n, abs(X), '.', 'MarkerSize', 20)
hold off
legend('Z', 'X')
title({'DFT(z) vs DFT(x)', 'centrada em Fa/2'});
xlabel('Amostras');
ylabel('Amplitude');

subplot(224)
stem(-L*N/2:L*N/2-1, abs(fftshift(Z)), 'x', 'MarkerSize', 10)
hold on
stem(-N/2:N/2-1, abs(fftshift(X)), '.', 'MarkerSize', 20)
hold off
legend('Z', 'X')
title({'DFT(z) vs DFT(x)', 'centrada em 0'});
xlabel('Amostras');
ylabel('Amplitude');

%% Alínea c
% Número de símbolos por amostra e respetivo vetor de indices
M = 2;
m = 0:N/M-1;

% Downsampling (M símbolos por amostra)
y = x(1:M:end);
Y = fft(y);

figure(6)
subplot(2,2,[1,2])
stem(n, x, 'x')
hold on
stem(m, y, 's')
hold off
title('Sinal y vs x')
xlabel('Indíce das amostras')
ylabel('Amplitude')

subplot(223)
stem(m, abs(Y), 'x', 'MarkerSize', 5)
hold on
stem(n, abs(X), '.', 'MarkerSize', 20)
hold off
legend('Y', 'X')
title({'DFT(y) vs DFT(x)', 'centrada em Fa/2'});
xlabel('Amostras');
ylabel('Amplitude');

subplot(224)
stem(-N/(M*2):N/(2*M)-1, abs(fftshift(Y)), 'x', 'MarkerSize', 5)
hold on
stem(-N/2:N/2-1, abs(fftshift(X)), '.', 'MarkerSize', 20)
hold off
legend('Y', 'X')
title({'DFT(y) vs DFT(x)', 'centrada em 0'});
xlabel('Amostras');
ylabel('Amplitude');

%% Exercício 5
% Frequência de amostragem da sinudoide anterior
Fa_x = 200;

%% Alínea a)
% Dados da sinusoide x[n]: N(número de amostras), k(valor da risca no espectro)
N = 100;
k = 40;

% Frequência [Hz] da sinusoide da sequência x[n]
f_x = Fa_x/N * k

%% Alínea b)
% Dados da sinusoide y[n]: N(número de amostras), k(valor da risca no espectro)
N = 100/M;
k = 10;

% Frequência de amostragem passa para metade porque o número de samples
% também Fa = Fs * Num_samples
Fa_y = Fa_x / M;

% Frequência [Hz] da sinusoide da sequência y[n]
f_y = Fa_y/N * k

%% Alínea c)
% Dados da sinusoide z[n]: N(número de amostras), k1,k2(valores da riscas no espectro)
N = 100*L;
k1 = 40;
k2 = 60;

% Frequência de amostragem passa para o dobro porque o número de samples
% também Fa = Fs * Num_samples
Fa_z = L * Fa_x;

% Frequência [Hz] da sinusoide da sequência z[n]
f_z1 = Fa_z/N * k1
f_z2 = Fa_z/N * k2
