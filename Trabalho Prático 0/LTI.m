%% Revisões: Sistemas Lineares Invariantes no Tempo (LTI)
clear; clc; close all;

%% a)
h.num = [1 2 1] * 0.443 * 10^-3;
h.den = [1 -1.94 0.94];

h.num_roots = roots(h.num);
h.den_roots = roots(h.den);

%% b)
if(isempty(h.den_roots))
    fprintf('O filtro é FIR\n')
else
    fprintf('O filtro é IIR\n')
end;

%% c)
[h.r, h.p, h.k] = residuez(h.num, h.den);

fprintf(['H\nr: [' strrep(sprintf(' %1.6f,',  h.r), ',', '' ) ']\n']);
fprintf(['p: [' strrep(sprintf(' %1.6f,',  h.p), ',', '' ) ']\n']);
fprintf(['k: [' strrep(sprintf(' %1.6f,',  h.k), ',', '' ) ']\n\n']);


if(find(round(h.den_roots, 4) == 1))
    fprintf('O sistema é marginalmente estável.\nPossui um pólo em 1\n');
end;

%% d)
zplane(h.num_roots, h.den_roots)

if(find(round(h.num_roots, 4) == -1))
    fprintf('Filtro passa baixo\n')
elseif(find(round(h.num_roots, 4) == 1))
    fprintf('Filtro passa alto\n')
end;
%% e)
x.num = [1 0 0 0 0 -1]
x.den = [1 -1]

[x.r, x.p, x.k] = residuez(x.num, x.den)

fprintf(['X\nr: [' strrep(sprintf(' %1.6f,',  x.r), ',', '' ) ']\n']);
fprintf(['p: [' strrep(sprintf(' %1.6f,',  x.p), ',', '' ) ']\n']);
fprintf(['k: [' strrep(sprintf(' %1.6f,',  x.k), ',', '' ) ']\n\n']);

y.num = conv(h.num, x.num);
y.den = conv(h.den, x.den);

[y.r, y.p, y.k] = residuez(y.num, y.den);

fprintf(['Y\nr: [' strrep(sprintf(' %1.6f,', y.r), ',', '' ) ']\n']);
fprintf(['p: [' strrep(sprintf(' %1.6f,',  y.p), ',', '' ) ']\n']);
fprintf(['k: [' strrep(sprintf(' %1.6f,',  y.k), ',', '' ) ']\n\n']);


%% 
d1 = 17.15;
d2 = 34.3;

v = 343;

delta_t1 = 2 * d1 / v
delta_t2 = 2 * d2 / v

Fa = 44100

N1 = Fa * delta_t1
N2 = Fa * delta_t2
