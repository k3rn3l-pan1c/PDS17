%% Revisões: Sistemas Lineares Invariantes no Tempo (LTI)
clear; clc; close all;

%% Ex1

% Sistema LTI
% Equação Diferença
% y[n] = 0.443 × 10^−3 x[n] + 0.886 × 10^−3 x[n − 1] + 0.443 × 10^−3 x[n − 2] 
%        + 1.94y[n − 1] − 0.94y[n − 2]

%% a)

% Da equação diferença retiramos a seguinte função de transferência:
h1.num = [1 2 1] * 0.443 * 10^-3;
h1.den = [1 -1.94 0.94];

% Raízes do numerador (zeros) e denominador (pólos)
h1.num_roots = roots(h1.num);
h1.den_roots = roots(h1.den);

%% b)

% O sistema é FIR se tiver todos os pólos na origem. Caso contrário é IIR
if(isempty(h1.den_roots))
    fprintf('O sistema é FIR\n')
else
    fprintf('O sistema é IIR\n')
end;

%% c)

% Separação em frações parciais na forma 
%    B(z)       r(1)               r(n)
%    ---- = ------------ +...  ------------ + k(1) + k(2)z^(-1) ...
%    A(z)   1-p(1)z^(-1)       1-p(n)z^(-1)
% 
% Caso existam pólos com multiplicidade diferente de 1, tal que se
% verifique P(j) = ... = P(j+m-1) existem termos na forma:
%    R(j)              R(j+1)                      R(j+m-1)
%   -------------- + ------------------   + ... + ------------------
%    1 - P(j)z^(-1)   (1 - P(j)z^(-1))^2           (1 - P(j)z^(-1))^m

[h1.r, h1.p, h1.k] = residuez(h1.num, h1.den);

 % Onde R é o resíduo, P os pólos e K os termos diretos
fprintf(['H1\nr: [' strrep(sprintf(' %1.6f,',  h1.r), ',', '' ) ']\n']);
fprintf(['p: [' strrep(sprintf(' %1.6f,',  h1.p), ',', '' ) ']\n']);
fprintf(['k: [' strrep(sprintf(' %1.6f,',  h1.k), ',', '' ) ']\n\n']);

if(find(round(h1.den_roots, 4) > 1))
    fprintf('O sistema é instável.\nPossui um pólo com amplitude superior a 1\n');
elseif(find(round(h1.den_roots, 4) == 1))
    fprintf('O sistema é marginalmente estável.\nPossui um pólo com amplitude 1\n');
else
    fprintf('O sistema é estável.\Todos os pólos possuem amplitude inferior a 1\n');
end;

%% d)

% Mapa de pólos e zeros
% x - pólo 
% o - zero
figure(1)
zplane(h1.num_roots, h1.den_roots)

% Considerando que o sistema só pode ser passa-baixo ou passa-alto, se
% tiver um zero em -1 é passa baixo porque atenua todas as frequências
% superiores a pi (Fa/2) porque introduz um zero no diagrama de Bode a essa
% frequência. Se tiver um zero em 1 (f=0), introduz um zero no diagrama de
% Bode a DC.
if(find(round(h1.num_roots, 4) == -1))
    fprintf('O sistema representa um filtro passa baixo\n')
elseif(find(round(h1.num_roots, 4) == 1))
    fprintf('O sistema representa um filtro passa alto\n')
end;

%% e)

% Entrada: x[n] = u[n] − u[n − 5] = 1/ (1 - z^-1) - z^-5 * 1/(1-z^-1)
%               = (1 - z^-5)/(1 - z^-1)

x.num = [1 0 0 0 0 -1];
x.den = [1 -1];

% A decomposição em frações parciais da entrada
[x.r, x.p, x.k] = residuez(x.num, x.den);

fprintf(['X\nr: [' strrep(sprintf(' %1.6f,',  x.r), ',', '' ) ']\n']);
fprintf(['p: [' strrep(sprintf(' %1.6f,',  x.p), ',', '' ) ']\n']);
fprintf(['k: [' strrep(sprintf(' %1.6f,',  x.k), ',', '' ) ']\n\n']);

% Convolucionando a entrada com a função de transferência (numerador com
% numerador e denominador com denominador), obtemos a saída
y.num = conv(h1.num, x.num);
y.den = conv(h1.den, x.den);

% Decomposição da saída em frações parciais
[y.r, y.p, y.k] = residuez(y.num, y.den);

fprintf(['Y\nr: [' strrep(sprintf(' %1.6f,', y.r), ',', '' ) ']\n']);
fprintf(['p: [' strrep(sprintf(' %1.6f,',  y.p), ',', '' ) ']\n']);
fprintf(['k: [' strrep(sprintf(' %1.6f,',  y.k), ',', '' ) ']\n\n']);

%% f) [ANALÌTICO]
aux.num = [1 2 1] * 0.443 * 10^-3;
aux.den = conv([1 -1.94 0.94], [1 -1]);

[aux.r, aux.p, aux.k] = residuez(aux.num, aux.den);

fprintf(['Aux\nr: [' strrep(sprintf(' %1.6f,', aux.r), ',', '' ) ']\n']);
fprintf(['p: [' strrep(sprintf(' %1.6f,',  aux.p), ',', '' ) ']\n']);
fprintf(['k: [' strrep(sprintf(' %1.6f,',  aux.k), ',', '' ) ']\n\n']);

% Error using residuez (line 60)
% First coefficient in A vector must be non-zero.
% [aux2.r, aux2.p, aux2.k] = residuez([0 1], conv([1 2 1], [0 -1]));
% 
% fprintf(['Aux2\nr: [' strrep(sprintf(' %1.6f,', aux2.r), ',', '' ) ']\n']);
% fprintf(['p: [' strrep(sprintf(' %1.6f,',  aux2.p), ',', '' ) ']\n']);
% fprintf(['k: [' strrep(sprintf(' %1.6f,',  aux2.k), ',', '' ) ']\n\n']);

%% Ex 2
h2.p = 0.5 * exp(1j*pi*[1 -1]/4);  % pólos
h2.z = [1 -1];                     % zeros
h2.k = 0.5;                        % ganho

%% a) [ANALÌTICO]

% H(z) = K * [ (1 - h2.z(1)*z^-1)/(1 - h2.p(1)*z^-1) + (1 - h2.z(2)*z^-1)/(1 - h2.p(2)*z^-1) ]

% Numerador e denominador da função de transferência
h2.num = conv([1 -h2.z(1)], [1 -h2.z(2)]);
h2.den = conv([1 -h2.p(1)], [1 -h2.p(2)]);

%% b)

% Mapa de pólos e zeros
figure(2)
zplane(h2.num, h2.den)
title('Poles and Zeros of H_2(z)');

% No domínio da frequência
figure(3)
freqz(h2.k * h2.num, h2.den)

%% c

% Função de transferência em fracções parciais (equivalente à resposta
% impulsional porque o sistema é LTI)
[h2.r, h2.p, h2.K] = residuez(h2.k * h2.num, h2.den)

fprintf(['H2\nr: [' strrep(sprintf(' %1.6f %+1.6fj,', real(h2.r), imag(h2.r)), ',', '' ) ']\n']);
fprintf(['p: [' strrep(sprintf(' %1.6f %+1.6fj,', real(h2.p), imag(h2.p)), ',', '' ) ']\n']);
fprintf(['k: [' strrep(sprintf(' %1.6f %+1.6fj,', real(h2.K), imag(h2.K)), ',', '' ) ']\n\n']);

% O resto é analítico

%% d
h2.step.num = 5;
h2.step.den = conv([1 -0.353553], [1 -1]);

% Em frações parciais
[h2.step.r, h2.step.p, h2.step.k] = residuez(h2.step.num, h2.step.den);

fprintf(['H2 step response\nr: [' strrep(sprintf(' %1.6f,', h2.step.r), ',', '' ) ']\n']);
fprintf(['p: [' strrep(sprintf(' %1.6f,',  h2.step.p), ',', '' ) ']\n']);
fprintf(['k: [' strrep(sprintf(' %1.6f,',  h2.step.k), ',', '' ) ']\n\n']);

%% e [DEMONSTRAÇÂO]

%% f

%% Exercício 3

% y[n] = x[n] + 0.8x[n − N_1 ] + 0.8^2 x[n − N_2 ]
%

%% a) [ANALÌTICO]

% Não sei a ordem dos coeficientes em z do polinómio. Como a função de
% transferência não apresenta recursividade com a saída, não tem pólos,
% logo é FIR

%% b)

% Distâncias do microfone ás paredes
d1 = 17.15;
d2 = 34.3;

% Velocidade de propagação do som no ar
v = 343;

% Frequência de amostragem
Fa = 44100;

% Considerando v = d/t e a distância de ida e de volta
delta_t1 = 2 * d1 / v;
delta_t2 = 2 * d2 / v;

% Cálculo do atraso dos ecos em número de amostras
N1 = Fa * delta_t1
N2 = Fa * delta_t2

%% c) [ANALÌTICO]

% Aplicando a transformada z ao sinal, colocamos o x(z) em evidência e em
% seguida o termo a*z^-D, que nos permite reescrever o resto da função como
% y.
% Após alguma matemática trivial, obremos y[n] = x[n] + a*y[n-D];

%% d) [ANALÌTICO]

%% e)

%% f)

%% Exercício 4

%% a)
h4.num = [ 1 2 2 1];
h4.den = 1;

% Raízes do numerador (zeros) e denominador (pólos)
h4.num_roots = roots(h4.num);
h4.den_roots = roots(h4.den);

figure(4)
zplane(h4.num_roots, h4.den_roots)

%% b) [ANALÌTICO]

%% c)
figure(5)
freqz(h4.num, h4.den)

%% d)
n = 0:3;
h5.num = (-1).^n .* h4.num;
h5.den = 1;

% Raízes do numerador (zeros) e denominador (pólos)
h5.num_roots = roots(h5.num);
h5.den_roots = roots(h5.den);

figure(6)
zplane(h5.num_roots, h5.den_roots)

%% e) 
figure(7)
freqz(h5.num, h5.den)

% O filtro passa a ter o comportamento dual: como era um passa baixo passou
% a ser um passa alto e vice-versa
