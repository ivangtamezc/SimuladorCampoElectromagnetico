%Entregable 3 - Campo Magnetico en el eje Z con RungeKutta de 4to orden

%Limpiamos el espacio de trabajo y la consola
clc
clear
%close all %Cierra todas las figuras abiertas

%Definicion Constantes
Mo =(4 * pi) * 1e-7;

%Inputs del usuario
MuiD = input("ingrese el momento magentico dipolar del iman: ");
MasaI = input("Ingrese la masa del iman: ");
Radio = input("Ingrese el radio del anillo: ");
I = input("Corriente electrica: ");

% MuiD = 70000;
% MasaI = 0.1;
% Radio = 0.8;
% I = 30.1;

%Creamos la malla de puntos
Vx = (-Radio:0.1:Radio);
Vy = (-Radio:0.1:Radio);
Vz = (0:0.1:3*Radio);
[X,Y,Z] = meshgrid(Vx, Vy, Vz);

%Creamos el vector del angulo theta
Theta = 0:0.01:2*pi;

%Parametros del iman
POZI = 2;%Posicion inicial del iman en el eje z
VOZ = 0;%Velocidad inicial del iman en el eje z

%Parametros del tiempo
Tf = 10;
Ti = 0;

%Creamos las posiciones del aro
PX = Radio * cos(Theta);
PY = Radio * sin(Theta);
PZ = zeros(size(PX));

%Parametros de Runge-Kutta de 4to orden
H = 0.001; %Tamaño del paso
T = Ti:H:Tf;

%Funcion que describe el cambio de campo
fAccel = @(T,X)((3 * MuiD * Mo * Radio^2 * I * X(1)) / (2 * MasaI * (Radio^2 + X(1)^2)^(5 / 3))) - 9.81;

% Condiciones iniciales
Posicion(1) = POZI;
Velocidad(1) = VOZ;
Aceleracion(1) = fAccel(T(1),Posicion(1));

% Metodo de Runge-Kutta de 4to orden
for i = 1:(length(T) - 1)
    %Calculamos lo coeficientes
    p1 = H * Velocidad(i);
    v1 = H * fAccel(T(i),Posicion(i));

    p2 = H * (Velocidad(i) + p1/2);
    v2 = H * fAccel(T(i) + H/2 ,Posicion(i) + H/2 * v1);

    p3 = H * (Velocidad(i) + p2/2);
    v3 = H * fAccel(T(i) + H/2 ,Posicion(i) + H/2 * v2) ;

    p4 = H * (Velocidad(i) + p3);
    v4 = H * fAccel(T(i) + H ,Posicion(i) + v3);

    Posicion(i+1) = Posicion(i) + (p1 + (2 * p2) + (2 * p3) + p4)/6;
    Velocidad(i+1) = Velocidad(i) + (v1 + (2 * v2) + (2 * v3) + v4)/6;
    Aceleracion(i+1) = fAccel(T(i),Posicion(i));

      if (Aceleracion(i+1) > -1e-6) && (Posicion(i+1) < 1.6e-6) && (Velocidad(i+1) > -1e-6)
          break;
      end

end

a = i + 1;

%Creamos la matriz del campo vacio
Bz = zeros(size(Z));

%Calculo del campo magnetico
for i = 1 :length(Vx)
    for j = 1 : length(Vy)
        for k = 1 : length(Vz)
            Bz(i,j,k) = (Mo * Radio^2 * I) / (2 * ( Radio^2 + Z(i,j,k)^2)^(3 / 2));
        end
    end
end

figure(Name="Campo Magnetico")
plot3(PX,PY,PZ)
hold on
quiver3(X,Y,Z,zeros(size(Bz)),zeros(size(Bz)),Bz,0.5)
axis equal
hold off

figure(Name="Graficas Runge-Kutta 4to Orden");
subplot(3, 1, 1);
plot(T(1:a), Posicion, 'b');
xlabel('Tiempo');
ylabel('Posición del imán');
title('Posición del imán en el eje z');

subplot(3, 1, 2);
plot(T(1:a), -Velocidad, 'r');
xlabel('Tiempo');
ylabel('Velocidad del imán');
title('Velocidad del imán en el eje z');

subplot(3, 1, 3);
plot(T(1:a), -Aceleracion, 'g');
xlabel('Tiempo');
ylabel('Aceleración del imán');
title('Aceleración del imán en el eje z');

%Creamos la figura 
figure(Name="Animacion del movimiento");

%Creamos subplots para ver diferentes angulos
subplot(1, 2, 1);

%Ploteamos el aro
Aro = plot3(PX,PY,PZ);
hold on

%Damos formato a los ejes
xlabel("Eje X")
ylabel("Eje Y")
zlabel("Eje Z")

%Limitamos los ejes
xlim([-Radio , Radio]);
ylim([-Radio , Radio]);
zlim([0 , max(Posicion)]);

%Ploteamos la posicicion de la gondola
gondola1 = plot3(0,0,Posicion(1),LineWidth=5,Color="b",Marker="o");

%Ajustamos el POV
view([180,180,20])
title('Vista en XYZ');

subplot(1, 2, 2);
%Ploteamos el aro
Aro = plot3(PX,PY,PZ);
hold on

%Damos formato a los ejes
xlabel("Eje X")
ylabel("Eje Y")
zlabel("Eje Z")

%Limitamos lo ejes
xlim([-Radio , Radio]);
ylim([-Radio , Radio]);
zlim([0 ,  max(Posicion)]);

%Ploteamos la posicicion de la gondola
gondola = plot3(0,0,Posicion(1),LineWidth=5,Color="b",Marker="o");

%Ajustamos el POV
view([180,0,0])
title('Vista en el Eje YZ');


%Actualizamos la posicicion de la gondola
for (i = 1 : length(T(1:a)))

    %Actualizamos la posicicion de la gondola
    set(gondola,"ZData",Posicion(i));
    set(gondola1,"ZData",Posicion(i));

    %Creamos un tiempo de espera para que se genere la animacion
    pause(0.001);

end