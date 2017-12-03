%% Exercício 1
clear;
clc;
close all;

%% a)
% Load the image and show
% Não posso usar im2double porque gera-me double \in [0, 1] e os pesos não
% estão nessa gama
image = double(imread('ImagensA/airfield02g.tif')) - 128;
figure(1)
imshow(image, [])
title('Original Image')

% Mascara
DCT_block = dctmtx(8);

% Função da DCT
my_dct = @(block_struct) DCT_block * block_struct.data * DCT_block';

% Aplicar a DCT bloco a bloco (8x8)
imDCT = blockproc(image, [8 8], my_dct);

% Imagem da DCT
figure(2)
imshow(imDCT, [])
figure(2)
title('DCT Transform')

%% b)
% Mascara de pesos a usar
steps = load('passos.txt');

% Função de DCT com os pesos
my_Q = @(block_struct) round(my_dct(block_struct) ./ steps);

% Quantificador
imQ = blockproc(image, [8 8], my_Q);

figure(3)
imshow(imQ, [])
title('DCT + Quantification')

%% c)
my_idct = @(block_struct) DCT_block' * block_struct.data * DCT_block;
my_Qinv = @(block_struct) my_idct(struct('data', block_struct.data .* steps));

% Desquantificador
imIDCT = blockproc(imQ, [8 8], my_idct);

figure(4)
imshow(imIDCT, [])
title('IDCT')

imQinv = blockproc(imQ, [8 8], my_Qinv);

figure(5)
imshow(imQinv, [])
figure(5)
title('IDCT + Q^{-1}')

%% d)
% Númeo de indices a considerar no cálculo da DCT
L = numel(steps);

% Empty vectors
ratio = zeros(1, L);
squareError = zeros(1, L);
Lvector = zeros(1, 8*8);

% Mostrar imagens?
verbose_fig = 0;

for k = 1:L
    Lvector(k) = 1;
    
    % Mascara zigzag para obter os pesos para cada pixel
    zigzagMask = izigzag(Lvector, 8, 8);
    
    % Quantificador em zigzag
    my_zigzagQ = @(block_struct) round(my_dct(block_struct) ./ steps .* zigzagMask );

    % Aplicar a Quantificação com zigzag bloco a bloco (8x8)
    imZigZagQ = blockproc(image, [8 8], my_zigzagQ);

    % Imagem quantizada com zigzag
    if(verbose_fig)
        figure(6)
        imshow(imZigZagQ, [])
        figure(6)
        title({'DCT + Quantization (with zigzag)', ['L= ' num2str(k)]})
    end;
    
    % Descodificador com desquantização em zigzag
    my_zigzagQinv = @(block_struct) my_idct(struct('data', block_struct.data .* steps .* zigzagMask));

    % Aplicar a desquantificação com zigzag bloco a bloco (8x8)
    imZigZagQinv =  blockproc(imZigZagQ, [8 8], my_zigzagQinv);

    % Imagem desquantizada com zigzag
    if(verbose_fig)
        figure(7)
        imshow(imZigZagQinv, [])
        figure(7)
        title({'IDCT + Q^{-1} (with Zigzag)', ['L= ' num2str(k)]})
    end;
    
    % O racio é o número de coeficientes úteis (diferentes de 0) e o número
    % de coeficientes que a mascara possui
    ratio(k) = (L - k) ./ L*100;
    
    % O erro médio quadrado é o somatório das diferenças ao quadrado da
    % imagem original e da imagem comprimida, a dividir pelo número de
    % elementos da imagem
    squareError(k) = sum(sum((image - imZigZagQinv).^2)) / numel(image);
end;

figure(8)
plot(ratio)
hold on
plot(squareError)
hold off
legend('Rácio de Compressão', 'Erro Médio Quadrático')
title('Rácio e Erro Médio Quadrático em função do Nº de coeficientes')
xlabel('Número de coeficientes a serem usados')
ylabel('Rácio / Erro Médio Quadrático')

%% e
% Númeo de indices a considerar no cálculo da DCT
L = numel(steps);

% Mostrar imagens?
verbose_fig = 0;

% Vetor de erros médios quadráticos
squareError = zeros(1, L);
imZigZagQinv = zeros(size(image));

for k = 1:L
    % Vector de indices 
    Lvector = zeros(1, 8*8);
    Lvector(k) = 1;
    
    % Mascara zigzag para obter os pesos para cada pixel
    zigzagMask = izigzag(Lvector, 8, 8);
    
    % Quantificador em zigzag
    my_zigzagQ = @(block_struct) round(my_dct(block_struct) ./ steps) .* zigzagMask ;

    % Aplicar a Quantificação com zigzag bloco a bloco (8x8)
    imZigZagQ = blockproc(image, [8 8], my_zigzagQ);

    % Imagem quantizada com zigzag
    if(verbose_fig)
        figure(9)
        imshow(imZigZagQ, [])
        figure(9)
        title({'DCT + Quantization (with zigzag)', ['L= ' num2str(k)]})
    end;
    
    % Descodificador com desquantização em zigzag
    my_zigzagQinv = @(block_struct) my_idct(struct('data', block_struct.data .* steps .* zigzagMask));

    % Aplicar a desquantificação com zigzag bloco a bloco (8x8)
    imZigZagQinv = imZigZagQinv + blockproc(imZigZagQ, [8 8], my_zigzagQinv);

    % Imagem desquantizada com zigzag
    if(verbose_fig)
        figure(10)
        imshow(imZigZagQinv, [])
        figure(10)
        title({'IDCT + Q^{-1} (with Zigzag)', ['L= ' num2str(k)]})
    end;
    
    % O racio é o número de coeficientes úteis (diferentes de 0) e o número
    % de coeficientes que a mascara possui
    ratio(k) = (L - k) ./ L*100;
    
    % O erro médio quadrado é o somatório das diferenças ao quadrado da
    % imagem original e da imagem comprimida, a dividir pelo número de
    % elementos da imagem
    squareError(k) = sum(sum((image - imZigZagQinv).^2)) / numel(image);
end;
    
figure(11)
plot(ratio)
hold on
plot(squareError)
hold off
legend('Rácio de Compressão', 'Erro Médio Quadrático')
title('Rácio e Erro Médio Quadrático em função do Nº de coeficientes')
xlabel('Número de coeficientes a serem usados')
ylabel('Rácio / Erro Médio Quadrático')

%% f) 
% Se tivesse de mostrar miniaturas de imagens, ao fazer codificação
% progressiva é possível obter rapidamente uma primeira versão de todas as
% imagens sem calcular a DCT completa. De seguida, calculo a DCT bloco a
% bloco para outro pixel para as imagens todas e sobreponho à imagem
% original. Assim consigo mostrar todas as imagens mais rapidamente e
% progressivamente aumento a qualidade delas. Este método possui vantagens
% quando precisa de ser aplicado a uma situação em que tem de gerar muitas
% miniaturas (uma miniatura não precisa que a imagem possua tanta
% qualidade, logo não precisa de usar tantos coefcientes no cálculo da DCT)
