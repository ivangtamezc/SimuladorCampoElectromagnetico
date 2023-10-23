%Entregable 1 - Reto

%Limpieza de la consola y espacio de trabajo
clc
clear
close all

%Definicion Constantes
Mo = 4 * pi * 1e-7;

%Inputs del usuario
dS = input("Ingrese el numero de elementos de carga: ");
Radio = input("Ingrese el radio del anillo: ");
I = input("Corriente electrica: ");
Malla = input ("Ingrese el tamaño de la malla: ");

%Definicion de la malla de vectores
Vx = linspace(-Malla,Malla,dS);
Vy = linspace(-Malla,Malla,dS);
Vz = linspace(-Malla,Malla,dS);

%Definicion de la malla
[X,Y,Z] = meshgrid(Vx,Vy,Vz);

%Ecuaciones
dtheta = 2*pi/dS;
Theta = 0:dtheta:(2*pi);

%Posicion del aro
PX = Radio*sin(Theta);
PY = Radio*cos(Theta);
PZ = zeros(size(PX));

%Diferenciales de posicion en el aro
DYi = X*dtheta;
DXi = -Y*dtheta;
DZi = zeros(size(DXi));

%Declaracion de los vectores vacios
Bx = zeros(size(X));
By = zeros(size(Y));
Bz = zeros(size(Z));

%Constante de la integral
Cte = Mo* I / (4 * pi);

for i=1:length(Vx)
    for j=1:length(Vy)
        for k=1:length(Vz)
                
            %Distancia en componentes de la malla al punto
            Dx(i,j,k) = X(i,j,k) - DXi(i,j,k);
            Dy(i,j,k) = Y(i,j,k) - DYi(i,j,k);
            Dz(i,j,k) = Z(i,j,k) - DZi(i,j,k);
            

            %Distancia total de la malla al punto
            Dt(i,j,k) = sqrt(Dx(i,j,k)^2 + Dy(i,j,k)^2 + Dz(i,j,k)^2);
  
            %Calculo del campo
            Bx(i,j,k) = Cte * (DYi(i,j,k)* Dz(i,j,k)/ Dt(i,j,k)^2) + Bx(i,j,k);
            By(i,j,k) = -Cte * (DXi(i,j,k)* Dz(i,j,k)/ Dt(i,j,k)^2)+ By(i,j,k); 
            Bz(i,j,k) = Cte * ((DXi(i,j,k)* Dy(i,j,k) - DYi(i,j,k)* Dx(i,j,k))/ Dt(i,j,k)^2) + Bz(i,j,k);

        end
    end
end

%Graficar
figure("Name", "Resultados Streamslice");
plot3(PX, PY, PZ , 'b', 'LineWidth', 2)
hold on
quiver3(X, Y, Z, Bx, By, Bz, 1.5,"filled");
title("Campo magnetico de un aro visto desde el plano XY");
xlabel("Eje X");
ylabel("Eje Y");
zlabel("Eje Z");