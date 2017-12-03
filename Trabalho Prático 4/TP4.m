%% Trabalho Prático 4 - Processamento Imagem: Filtros
clear;
clc;
close all;

% (Dis)Enable some pictures print 
verbose_pic = 0;

%% Ex 1
% Filtro Gaussiano
L = 3:2:15;
sigma = [ 0.5 0.8 1 1.2 1.5 1.8 2 ];

for k = 1:length(L);
    % Escolher k-th elementos de ambos os arrays
    H1 = fspecial('gaussian', L(k), sigma(k));

    figure(1)
    suptitle('Gaussian Filter dependence of \sigma')
    subplot(floor(sqrt(length(L)))+1, floor(sqrt(length(L)))+1, k)
    freqz2(H1)
    axis([-1 1 -1 1 0 1])
    title(['\sigma = ' num2str(sigma(k)) ])
end;

% Quanto maior o sigma, maior o desvio padrão, ou seja, a variação entre
% frequências adjacentes é mais abrupta. Isto causa com que o filtro tenha
% maior magnitude perto da sua média (0, neste caso) e rapidamente
% decresça atingido magnitudes muito baixas. Analisando o caso de sigma =
% 0.5, verificamos que ainda existe amplitude não nula a Fa/2 mas para
% sigma = 2 em Fa/4 já possuimos amplitude nula. 
% Na prática isto significa que aumentar o sigma torna o filtro mais
% selectivo às baixas frequências (melhora o seu fator de qualidade).
% Este filtro tem um comportamento passa baixo e quanto maior o sigma, mais
% amaciada fica a imagem porque o filtro suaviza as zonas de transições
% abruptas, sendo removidas as zonas de alto contraste, ou seja, as baixas
% frequências são preservadas enquanto as altas frequências são eliminadas.
% A dependência do fator sigma no filtro é contrário ao que se verifica
% para o tempo. No tempo se o sigma for mais pequeno, o filtro é mais curto
% no tempo 

%% Ex2
% Filtro Realce (sharpen)
h = [-0.5; 2; -0.5];
H2 = h*h';
H2(2,2) = 3;
figure(2)
freqz2(H2)
axis([-1 1 -1 1 0 2])
title('Filtro de Realce')

% O filtro atenua fortemente as baixas frequências e amplifica as altas
% frequências. Assim é bom para detetar arestas numa imagem porque sempre
% que existam transições de cor/intensidade o filtro amplifica essas
% transições. Se as transições entre pixeis forem suaves, ou seja, os
% pixeis tenham aproximadamente a mesma cor, a frequência da imagem nesses
% pixeis á baixa, pelo que essas frequências/diferenças são eliminadas pelo
% filtro. Desse modo, diz-se que o filtro é um filtro de Realce/detetor 
% de arestas porque realça as arestas das formas na imagem, que são
% componenetes de alta frequência. As arestas são detetadas quer segundo o 
% eixo vertical quer segundo o eixo horizontal, visto ser simétrico em
% relação a ambos os eixos. 

%% Ex3
% Read images
files = dir('peppers/pepn*.tif');
files(end+1) = dir('peppers/pepb*.tif');

% Read and show original image
image.peppers = im2double(imread('peppers/peppers.tif'));

figure(3)
subplot(sqrt(length(files)+1), sqrt(length(files)+1), 1)
imshow(image.peppers)
title('Peepers Original')

MSE = zeros(1, length(files));
for k = 1 : length(files)
    % Get images with noise and convert to double
    noise_peppers = im2double(imread(['peppers/' files(k).name]));
    
    % show image
    subplot(sqrt(length(files)+1), sqrt(length(files)+1), k+1)
    imshow(noise_peppers, [])
    title(files(k).name);
    
    % Calculate the Mean Square Error
    MSE(k) = sum(sum(image.peppers - noise_peppers).^2)/numel(image.peppers);
end;

figure(4)
plot(MSE, '-x')
set(gca,'XTick',1:3)
set(gca,'XTickLabel',{files.name})
title('Mean Square Error')
% xlabel('Images')
ylabel('MSE')

%% Ex4
image.pepnoise1 = double(imread('peppers/pepnoise1.tif'));
figure(5)
subplot(floor(sqrt(length(L))), ceil(sqrt(length(L)))+1, 1)
imshow(image.pepnoise1, [])
title('Pepnoise1 Original')
colorbar

% Filtrar coms os vários coeficientes do filtro Gaussiano
for k = 1:length(L);
    % Escolher k-th elementos de ambos os arrays
    H1 = fspecial('gaussian', L(k), sigma(k));

    im = imfilter(image.pepnoise1, H1, 'conv', 'same');
    % show image
    figure(5)
    subplot(floor(sqrt(length(L))), ceil(sqrt(length(L)))+1, k+1)
    imshow(im, [])
    colorbar
    title(['Guassian Filter applied to Pepnoise | \sigma=' num2str(sigma(k))]);
    
    % Calculate the Mean Square Error
    MSE(k) = sum(sum(im - image.pepnoise1).^2)/numel(im);
end;


figure(6)
plot(MSE, '-x')
set(gca,'XTick', 1:length(L))
set(gca,'XTickLabel',num2str(sigma'))
title('Mean Square Error')
xlabel('\sigma')
ylabel('MSE')

% A imagem com mais qualidade é a imagem que é filtrada com sigma = 1.2
% No entanto, a imagem que possui um menor MSE é a imagem filtrada com um
% sigma = 2. Assim podemos concluir que não existe uma relação direta entre
% o MSE e a qualidade aparente da imagem que é percetida pelo olho humano,
% uma vez que as métricas não são coerentes


%% Ex5
% Ler imagem
image.pepblur = imread('peppers/pepblur.tif');

% Filtrar usando o H2 (com convulação)
image.pepblur_filtered = imfilter(double(image.pepblur), H2, 'conv', 'same');

figure(7)
imshow(image.pepblur_filtered, [])
figure(7)
title('Pepblur filtered with  a Sharpen FIlter')
colorbar

%% a) \alpha = 1

% Considerando o alfa = 1
alfa = 1;
out = alfa * image.pepblur_filtered + double(image.pepblur);

figure(8)
subplot(131)
imshow(image.pepblur, [])
title('Pepblur Original')
colorbar

subplot(132)
imshow(out, [])
title('Pepblur filtrado e quantizado para Double')
colorbar

subplot(133)
imshow(uint8(out))
title('Pepblur filtrado e quantizado para Uint8')
colorbar

% Comparando  os resultados das imagens podemos concluir que obtemos
% resultados melhores quando após a filtragem efetuamos a quantização para
% apenas 255 níveis, ao invés de ter todos os níveis como na imagem
% double

%% b) 
alfa = 0.1:0.1:2;
verbose_pic = 1
for k = 1 : length(alfa)
    out = alfa(k) * image.pepblur_filtered + double(image.pepblur);
    
    % Show pictures if desired
    if verbose_pic
        figure(8 + k)
        subplot(121)
        imshow(out, [])
        colorbar
        subplot(122)
        imshow(uint8(out))
        colorbar
        suptitle(num2str(alfa(k)))
    end;

    % Calculate the Mean Square Error
    MSE(k) = sum(sum(out - double(image.pepblur)).^2)/numel(image.pepblur);
end;

% Mean Square Error
figure(9 + verbose_pic *  k)
plot(MSE, '-x')
set(gca,'XTick', 1:length(alfa))
set(gca,'XTickLabel',num2str(alfa'))
title('Mean Square Error')
xlabel('\sigma')
ylabel('MSE')

% Este procedimento tem como objetivo tentar compensar o blur da imagem ao
% sobrepor uma percentagem da imagem original filtrada por um filtro
% detetor de arestas com a imagem blured. Assim é possível reduzir o
% efeitodo blur ao "acrescentar" contraste com as arestas dos objetos que
% forma removidas quando se efetuou o blur

%% Ex6
MSE_pepnoise = zeros(1, 3);

% Obter imagem original
image.pepnoise2.original = imread('peppers/pepnoise2.tif');

% Filtrar com H1 e calcular respetivo MSE
image.pepnoise2.H1 = imfilter(double(image.pepnoise2.original), H1, 'conv', 'same');
MSE_pepnoise(1) = sum(sum(double(image.pepnoise2.original) - image.pepnoise2.H1).^2)...
                  /numel(image.pepnoise2.original);

% Filtrar com H2 e calcular respetivo MSE
image.pepnoise2.H2 = imfilter(double(image.pepnoise2.original), H2, 'conv', 'same');
MSE_pepnoise(2) = sum(sum(double(image.pepnoise2.original) - image.pepnoise2.H2).^2)...
                  /numel(image.pepnoise2.original);

image.pepnoise2.medfilter = medfilt2(double(image.pepnoise2.original));
MSE_pepnoise(3) = sum(sum(double(image.pepnoise2.original) - image.pepnoise2.medfilter).^2)...
                  /numel(image.pepnoise2.original);

% Display images
figure(10 + verbose_pic *  k)
subplot(221)
imshow(image.pepnoise2.original, [])
title('Original')

subplot(222)
imshow(image.pepnoise2.H1, [])
title('H1')

subplot(223)
imshow(image.pepnoise2.H2, [])
title('H2')

subplot(224)
imshow(image.pepnoise2.medfilter, [])
title('Median Filter')

% Mean Square Error
figure(11 + verbose_pic *  k)
plot(MSE_pepnoise, '-x')
set(gca,'XTick', 1:3)
set(gca,'XTickLabel',{'H1', 'H2','Non-Linear Filtering'} )
title('Mean Square Error')
xlabel('\sigma')
ylabel('MSE')

disp('Ex6')
T = table(MSE_pepnoise', 'RowNames', {'H1', 'H2', 'MedianFilter'}, ...
          'Variablenames', {'MSE'})

%% 7
%% a) 
% out = filter(h, 1, filter(v, 1, im)')'
% 
% O comando acima implementa a decomposição de filtros para imagem 2D
% segundo as suas duas dimensões, criando filtros cascateados (1Dx1D)

%% b
% Filtro separável
v = ones(1, 3)' / 3;    % Filtro segundo as linhas
h = ones(1, 3)' / 3;    % FIltro segundo as colunas

% Filtro cascateado
image.peppers_separated_filt = filter(h, 1, filter(v, 1, image.peppers)')';

figure(12 + verbose_pic *  k)
imshow(image.peppers_separated_filt)
title('Peppers com Filtro Cascateável')

%% c
% Usando a single value decompostion obter as matrizes U,S e V para o
% filtro H1
[U, S, V] = svd(H1);

% O filtro é separável se a decomposição for reversível, ou seja, o produto
% destes vcetores originar a matrix do filtro incial
if (isequal(round(U*S*V', 5),H1))
    disp('O filtro H1 é separável')
else
    disp('O filtro H1 não é separável')
end;

% Usando a single value decompostion obter as matrizes U,S e V para o
% filtro H2
[U, S, V] = svd(H2);

% O filtro é separável se a decomposição for reversível, ou seja, o produto
% destes vcetores originar a matrix do filtro incial
if (isequal(round(U*S*V', 5),H2))
    disp('O filtro H1 é separável')
else
    disp('O filtro H2 não é separável')
end;
