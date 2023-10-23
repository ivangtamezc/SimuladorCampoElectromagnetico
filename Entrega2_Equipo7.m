%Entregable 2 - Reto

%-Limpieza de la consola y espacio de trabajo-
clc
clear
close all

%---------------Definicion Constantes---------
Mo =(4 * pi) * 1e-7;

%----------------Inputs del usuario-----------
dS = input("Ingrese el numero de elementos de carga: ");
Radio = input("Ingrese el radio del anillo: ");
I = input("Corriente electrica: ");
Malla = input ("Ingrese el tamaño de la malla: ");
nE = input ("Ingrese el número de espiras: ");
dE = input ("Ingrese la distancia entre los aros: ");

%------Definicion de la malla de vectores-----
Vx = linspace(-Malla,Malla,dS);
Vy = linspace(-Malla,Malla,dS);
Vz = linspace(-Malla,Malla,dS);

%-------------Definicion de la malla----------
[X,Y,Z] = meshgrid(Vx,Vy,Vz);

%----------------------Aro--------------------
aux = 1;
for q = 1:nE
    %---------------Obtener el ángulo-------------
    dtheta = (2*pi)/dS;
    Theta = 0:dtheta:(2*pi)-dtheta;
    %---------Obtener la posición del aro---------
    PX(aux:aux+dS-1) = Radio*cos(Theta);
    PY(aux:aux+dS-1) = Radio*sin(Theta);
    PZ(aux:aux+dS-1) = (-nE*dE/2)+dE*(q-1);
    %-----Diferenciales de posicion en el aro-----
    DXi(aux:aux+dS-1) = -PY(aux:aux+dS-1) * dtheta;
    DYi(aux:aux+dS-1) = PX(aux:aux+dS-1) * dtheta;
    aux = aux + dS;
end
DZi(1:nE*dS) = 0;

%---------------------Graficar---------------
figure("Name", "Bombilla");
quiver3(PX, PY, PZ, DXi, DYi, DZi, 0.5);
axis equal

%-----Declaracion de los vectores vacios-----
Bx = zeros(size(X));
By = zeros(size(Y));
Bz = zeros(size(Z));
Dx(1:length(Vx),1:length(Vy),1:length(Vz)) = 0;
Dy(1:length(Vx),1:length(Vy),1:length(Vz)) = 0;
Dz(1:length(Vx),1:length(Vy),1:length(Vz)) = 0;
Dt(1:length(Vx),1:length(Vy),1:length(Vz)) = 0;

%----------Constante de la integral----------
CteB = (Mo * I) / (4 * pi);

for i=1:length(Vz)
    for j=1:length(Vy)
        for k=1:length(Vx)
            for p=1:(nE*dS)
            %Distancia en componentes de la malla al punto
            Dx(i,j,k) = X(i,j,k) - PX(p);
            Dy(i,j,k) = Y(i,j,k) - PY(p);
            Dz(i,j,k) = Z(i,j,k) - PZ(p);

            %Distancia total de la malla al punto
            Dt(i,j,k) = sqrt(Dx(i,j,k).^2 + Dy(i,j,k).^2 + Dz(i,j,k).^2 + (0.2).^2);
  
            %Calculo del campo
            Bx(i,j,k) = (CteB * (DYi(p) * Dz(i,j,k)/ Dt(i,j,k).^3)) + Bx(i,j,k);
            By(i,j,k) = -(CteB * (DXi(p) * Dz(i,j,k)/ Dt(i,j,k).^3)) + By(i,j,k); 
            Bz(i,j,k) = (CteB * ((DXi(p) * Dy(i,j,k) - DYi(p)* Dx(i,j,k))/ Dt(i,j,k)^3)) + Bz(i,j,k);
            end
        end
    end
end

B = sqrt(Bx.^2 + By.^2 + Bz.^2);

%------------------Graficar-----------------
figure("Name", "Resultados quiver 3D");
plot3(PX, PY, PZ , 'b', 'LineWidth', 2)
hold on
quiver3(X, Y, Z, Bx, By, Bz,"filled");
title("Campo magnetico de un aro");
xlabel("Eje X");
ylabel("Eje Y");
zlabel("Eje Z");

%Graficar 2D YZ
origen=round(length(Vx)/2);
By_xz=squeeze(By(:,origen,:));
Bz_xz=squeeze(Bz(:,origen,:));
figure("Name", "Resultados Streamslice 2D del plano YZ");
plot(PY,PZ)
streamslice(Vy,Vz,By_xz',Bz_xz',2)
xlabel("Eje Z");
ylabel("Eje Y");