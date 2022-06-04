%% Proyecto -> Goniometro empleando la trasformada de Hough

% Grado y grupo: 7A
% Equipo: 2
    % Flores Esparza, María Yessenia
    % Herrera Avitia, María Fernanda
    % Laris Santos, Andrea
    % Serna de la Torre, Verónica


clc
clearvars

%% Cargar una imagen
I = rgb2gray(imread('Images/B5.jpg'));

% Filtro gaussiano, suavizar imagen 
B = imgaussfilt(I,1); 
    % Mandar imagen y sigma (desviacion esstándar, la cual al ser pequueña
    % elimina frecuencias altas

% Descartar intensidades mayores a 20
B = B<20;   
    % Considerar solo los puntos negros de la cinta empleada

% Binarizar imagen
BW = edge(B,'canny',[0.01,0.35]);
    % Detectar los bordes de la imagen, el rango depende de la imagen

% Gráficar imagen original, filtrada y binaria
figure(1);
subplot(1,3,1);imagesc(I); colorbar();title('Original Image'); colormap();
subplot(1,3,2);imagesc(B); colorbar();title('Smoothed Image');
subplot(1,3,3);imagesc(BW); colorbar();title('Binary Edges');

% % Código auxiliar para apreciar imagenes por separado
% figure(11);
% imagesc(BW); colorbar();title('');

%% Declarar arreglos a utilizar

[x, y] = size(I); % Dimensiones de la imagen 

% Vector para guardar valores de rho
d = round(sqrt(x^2 + y^2)); % Sacar hipotenusa de la imagen
                            % considerando distancia de una esquina a otra
rho = (-d:5:d);

% Vector para guardar valores de theta
theta = (-89:1:90);

% Hough (matriz de votos, [rho,theta])
H = zeros(length(rho),length(theta)); 

%% Llenar matriz de votaciones

for l= 1:x % renglones de la imagen
    for m= 1:y  % columnas de la imagen
        if(BW(l,m)==1) % Calcular solo para puntos blancos (bordes)
            for k =1:length(theta)%recorrimos todo theta para calcular rho
                r= round(l*cos(theta(k)*pi/180) + m*sin(theta(k)*pi/180));
                
                b=find(rho==r); 
                % Buscar en el vector rho, el valor obtenido para r y
                % obtenemos esa posicion
                
                H(b,k)=H(b,k)+1;
                % Acumular votos según los valores de theta y rho
               
            end
        end
    end
end

% Graficar la matriz votada -> Hoough (espacio paramétrico theta, rho)
figure(2)
imagesc(theta,rho,H);colorbar();title('Transforamda de Hough');
xlabel('\theta'),ylabel('\rho'); 

%% Encontrar picos
    % Primero se busca la casilla más votada por columnas (theta) con lo
    % cual se selecciona la linea más larga de todas aquellas con el mismo
    % ángulo. Despues, se encuentran los 2 valores de theta más votados. 
    % Al final se guardan las posiciones

mayor=0; % Variable auxiliar
Htest=H; % Matriz de apoyo para encontrar valores máximos para theta
posY=zeros(1,2); % Guardar las posiciones en theta para 2 picos
posX=zeros(1,2); % Guardar las posiciones en rho para 2 picos

for m=1:length(posY) % Rutina ejecutada 2 veces (2 números de picos)
    for k=1:size (H,1) % Recorrer renglones de Hough
        for l=1:size(H,2) % Recorrer columnas de Hough
            if(Htest(k,l)>mayor) % Encontrar el pico mayor
                mayor=Htest(k,l);
                posY(1,m) = k; % Guardar la posicion en y del pico
                posX(1,m) = l; % Guardar la posicion en x del pico
            end
        end
    end
    Htest(:,posX(1,m)-5:posX(1,m)+1)=0; % Llenar de ceros la columna (theta) donde se encontró
                          % el primer pico
    mayor=0;
end

% Graficar Hough señalando los picos encontrados 
% (las dos lineas más largas encontradas)
figure(3)
imshow(H,[],'XData',theta,'YData',rho,'InitialMagnification','fit');
title('Transforamda de Hough con picos identificados');colorbar();
xlabel('\theta'),ylabel('\rho'); 
axis on, axis normal, hold on;

picT = theta(posX(1,:)); picR = rho(posY(1,:));
plot(picT,picR,'s','color','red');

% picT, almacena los valores de theta de los picos encontrados
% picR, almacena los valores de rho de los picos encontrados

if(picT(:)==90)
    %% Encontrar lineas para "y" de la ecuación
        % Se despejó "y" para evitar ideterminaciones debido al ángulo
        % en caso de angulos de 90 grados

    xCar = 1:size(I,2); % Generar valores para x de 1 al ancho de la imagen
    yCar1 = zeros(1,length(xCar)); % Guardar valores de "x" para el primer pico
    yCar2 = zeros(1,length(xCar)); % Guardar valores de "x" para el segundo pico

    % Primer pico
    for m=1:length(xCar)
        yCar1(m)= (picR(1) - (xCar(m)*cos(picT(1)*pi/180))) / sin(picT(1)*pi/180);
    end

    % Segundo pico
    for m=1:length(xCar)
        yCar2(m)= (picR(1) - (xCar(m)*cos(picT(1)*pi/180))) / sin(picT(1)*pi/180);
    end

    xLine2=xCar'; % Transponer x para poder dibujar con line

    % Guardar en una matriz vertical los valores en "y" del primer y segundo pico
    yLine2(:,1) = yCar1(:);
    yLine2(:,2) = yCar2(:);

    figure(8)
    imagesc(I);colorbar();title('Lineas identificadas')
    axis on, axis normal, hold on;
    line(yLine2,xLine2,'LineWidth',5); 
    hold off;

else
    %% Encontrar lineas despejando "x" de la ecuación
        % Se despejó "x" para evitar ideterminaciones debido al ángulo

    yCar = 1:size(I,1); % Generar valores para y de 1 a la altura de la imagen
    xCar1 = zeros(1,length(yCar)); % Guardar valores de "x" para el primer pico
    xCar2 = zeros(1,length(yCar)); % Guardar valores de "x" para el segundo pico

    % Primer pico
    for m=1:length(yCar)
        xCar1(m)= (picR(1) - (yCar(m)*(sin(picT(1)*pi/180))))/cos(picT(1)*pi/180);
    end

    % Segundo pico
    for m=1:length(yCar)
        xCar2(m)= (picR(2) - (yCar(m)*(sin(picT(2)*pi/180))))/cos(picT(2)*pi/180);
    end

    yLine=yCar'; % Transponer y para poder dibujar con line

    % Guardar en una matriz vertical los valores en x del primer y segundo pico
    xLine(:,1) = xCar1(:); 
    xLine(:,2) = xCar2(:);

    % Graficar imagen original sobreponiendo las lineas encontradas
        % Para solucionar la referencia se mandan a la función line 
        % primero los valores en "x" y depues los valores en "y"
    figure(4)
    imagesc(I);colorbar();title('Lineas identificadas') % este mero
    axis on, axis normal, hold on;
    line(yLine,xLine,'LineWidth',5); 
    hold off;

    % % Código auxiliar para apreciar lineas con referencia según MatLab
    % figure(41)
    % imagesc(I);colorbar();title('Lineas identificadas sin girar')
    % axis on, axis normal, hold on;
    % line(xLine,yLine,'LineWidth',5); 
    % hold off;
end

%% Encontrar el ángulo entre lineas identificadas

sumA = sum(picT);
if(sumA>=0) 
    angulo = 180-sum(abs(picT));
else
    angulo = max(abs(picT))+ max(picT);
end

disp(['El ángulo encontrado es de ',num2str(angulo),' grados']);

title(['El ángulo encontrado es de ',num2str(angulo),' grados']);

