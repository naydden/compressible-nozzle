clc;
clear all;
tic
%% Entrada de datos. 
ra=0;
d=10^-6;
%%Gas : aire
cp=@(T) 1022-0.1626*T+3.5025*10^-4*T^2;
lambda=@(T) 3.807*10^-3+7.4*10^-5*T;
nu=@(T) (2.53928*10^-5*sqrt(T/273.15))/(1+(122/T));

%%Caso de tobera convergente-divergente
er=0.0001;
L=0.5;
N=555;
dz=L/N;
conv_div=1;
for i=1:(N+1)
    Rs(i)=1/(100*i*dz+1.5)+i*dz+0.015;
    D(i)=0.01;
    S(i)=pi*Rs(i)^2;
end
[MIN , Nmin]=min(Rs(:));
for i=1:(N)
    theta(i)=atan((Rs(i+1)-Rs(i))/dz);
    Al(i)=2*pi*Rs(i)*dz/cos(theta(i));
end
    theta(N+1)=atan((Rs(N+1)-Rs(N))/dz);
    Al(N+1)=2*pi*Rs(N+1)*dz/cos(theta(N+1));

Pin=5*100000;
Tin=600; 
R=287;
gamma=(cp(Tin))/(cp(Tin)-R);
Tw=300; %in K and cte

%% Obtencion de Datos a diferentes M de entarda

for Min=[0.01 0.03 0.04 0.05571];
Vin=Min*sqrt((gamma*R*Tin)); %m/s
ra=ra+1;
% Caso de tobera no-crítica (obtenido por pruebas)
                % Onda de choque se ve a producir despues del tubo
                % (no hay onda de choque)
                choq=N+5;
                %%Se llama al programa que realiza todos los calculos
                compresible;
                % Se guardan los resultados
                Pnc(ra,:)=P; Tnc(ra,:)=T; Vnc(ra,:)=V; Mnc(ra,:)=M;
                entrnc(ra,:)=entr; Sgennc(ra,:)=Sgen;
% Caso de tobera crítica
          if Min==0.05571
                %puntos donde se va a calcular la onda de choque
                %(Nmin=garganta)
                rw=1;
                w=round(linspace(Nmin+100,N,4));
                for choq=w
                %%Se llama al codigo que realiza todos los calculos para cada onda de
                %%choque
                compresible;
                Pc(rw,:)=P; Tc(rw,:)=T; Vc(rw,:)=V; Mc(rw,:)=M;
                entrc(rw,:)=entr; Sgenc(rw,:)=Sgen;
                rw=rw+1;
                end
          end
end
DataFINAL = struct('Pnc',Pnc,'Pc',Pc,'Tnc',Tnc,'Tc',Tc,'Vnc',Vnc,'Vc',Vc,'Mnc',Mnc,'Mc',Mc,'entrnc',entrnc,'entrc',entrc,'Sgennc',Sgennc,'Sgenc',Sgenc); 
%% Tratamiento de Datos
graficas(DataFINAL,L,N);
toc