%Clasificación no lineal de patrones 
%Juan Camilo Amaya
%Manuela Cerón Viveros
%Juan Daniel Muñoz

clear all
clc;
close all;
display('Punto 2 - Parcial 2');
Ne=20; %neuronas en capa oculta
N=2000; %Muestras 
in=2; %Neuronas de entrada
out=1; %Neuronas de salida
f=0.02; %Factor para función Tanh

r1=10; %Radio 1
r2=6; %Radio 2
PX1=ones(3,N); %Inicialización media luna positiva y target
PX2=-1*ones(3,N); %Inicialización media luna negativa y target
k1=0;
k2=0;
d=-4; %Distancia vertical entre medias lunas
r=8; %Distancia horizontal entre medias lunas

%Generación de Media luna positiva
while k1 <= N-1
    %Coordenadas aleatorias
    px1=randi([-100*r1,100*r1])/100;
    py1=randi([0,100*r1])/100;
    %Magnitud
    a=sqrt(px1^2+py1^2);
    %Verificación de rango
    if a<r1 && a>r2
        PX1(1,k1+1)=px1;
        PX1(2,k1+1)=py1;
        k1=k1+1;
  
    end
    
end
%Generación de Media luna negativa
while k2 <= N-1
    %Coordenadas
    px2=randi([-100*r1,100*r1])/100;
    py2=(randi([-100*r1,0])/100);
    %Magnitud
    b=sqrt(px2^2+py2^2);
    %Verificación de rango
    if b<r1 && b>r2
        PX2(1,k2+1)=px2+r;
        PX2(2,k2+1)=py2-d;
        k2=k2+1;
    end

end

P=[PX1 PX2]; %Vetor de medias lunas y target

figure(1)
hold on
plot(PX1(1,:),PX1(2,:),'.r')
hold on
plot(PX2(1,:),PX2(2,:),'.b')
title('Entrada para entrenamiento')

%%
%Entrenamiento

Wji=-0.5+(0.5).*rand(Ne,in);%matriz de pesos sinápticos entrada - capa oculta
Wkj=-0.5+(0.5).*rand(out,Ne);%matriz de pesos sinápticos capa oculta - capa salida 

PX=P;
for n=1:100;
    alphax=-0.0101*n+1.0101; %Variable entre 1- 1x10^-4
for k=1:N;
    Uj=Wji*PX(1:2,k);%entrada neta en los nodos ocultos
    Yj=tanh(f*Uj);%función de transferencia nodos ocultos
    Dy=f*(1-Yj.^2);
    Vk=Wkj*Yj;%entrada neta en los nodos de salida
    Zk=tanh(f*Vk);%función de transferencia nodos de salida
    Dz=f*(1-Zk.^2);
    e(:,k)=PX(3,k)-Zk; %error
    deltaS=e(:,k).*Dz;%cálculo del delta específico de salida
    
    for j=1:Ne;
        deltaO(j)=Dy(j).*deltaS'*Wkj(:,j);%cálculo del delta específico oculto
    end
    clear j
    
    for j=1:out; %Para una neurona de salida
        for i=1:Ne;
            Wkj(j,i)=Wkj(j,i)+alphax*deltaS(j).*Yj(i);%actualización de los pesos oculta-salida
        end
        clear i
    end
    clear j
    for j=1:Ne;%actualización de los pesos entrada-oculta
        for i=1:in;
            Wji(j,i)=Wji(j,i)+alphax*deltaO(j)*PX(i,k);
        end
        clear i
    end
    clear j
    
    er(n,k)=0.5*sum(e(:,k).^2); 
    clear Yj
    clear Zk
end
clear e
clear Uj
clear Vk
PX=P(:,randperm(N)); %Permutación de la entrada y target
end
erx=(sum(er,2))/N;%error cuadrático medio
figure(2)
plot(erx)
title('Error cuadrático Medio')
%%
%Validación individual
%MEdia luna positiva
maxX1=max(PX1(1,:));
minX1=min(PX1(1,:));
maxY1=max(PX1(2,:));
minY1=min(PX1(2,:));

%MEdia luna negativa
maxX2=max(PX2(1,:));
minX2=min(PX2(1,:));
maxY2=max(PX2(2,:));
minY2=min(PX2(2,:));

x=input('Coordenada en x=');
y=input('Coordenada en y=');

if (x<=minX1 || x>=maxX1) && (x<=minX2 || x>=maxX2) 
    l=1;
else
    l=0;
    
end

 if (y<=minY1 || y>=maxY1) && (y<=minY2 || y>=maxY2) 
     z=1;
 else
     z=0;
     
 end
 i=l+z;
 if i==0
     ve=[x;y];

Netaop=Wji*ve;
Yj=tanh(f*Netaop);
Netasp=Wkj*Yj;
Zk=round(tanh(f*Netasp))
 else
 disp('Entrada no válida');
   
 end
 
 
%%
%Validación de muestra completa

Nv=N*2;
PX1v=ones(3,Nv);
PX2v=-1*ones(3,Nv);
dv=-4;
k3=0;
k4=0;
%Generación de entrada de validación
while k3 <= Nv-1
    px1v=randi([-100*r1,100*r1])/100;
    py1v=randi([0,100*r1])/100;
    
    av=sqrt(px1v^2+py1v^2);
    
    if av<r1 && av>r2
        PX1v(1,k3+1)=px1v;
        PX1v(2,k3+1)=py1v;
        k3=k3+1;
    end
end
%  clear k1

while k4 <= Nv-1
    px2v=randi([-100*r1,100*r1])/100;
    py2v=(randi([-100*r1,0])/100);
    
    bv=sqrt(px2v^2+py2v^2);
    
    if bv<r1 && bv>r2
        PX2v(1,k4+1)=px2v+r;
        PX2v(2,k4+1)=py2v-dv;
        k4=k4+1;
    end
xlim([-20 20])
ylim([-20 20])
end
%  clear k2


Pv=[PX1v PX2v]; %Vector de entrada
PVAL=Pv(:,randperm(Nv)); %Permutación de la entrada para tener coordenadas aleatorias

h=0;
clear y;clear z;
for n=1:100;
    Ujv=Wji*PVAL(1:2,n);%entrada neta en los nodos ocultos
    Yjv=tanh(f*Ujv);%función de transferencia nodos ocultos
    Vkv=Wkj*Yjv;%entrada neta en los nodos de salida
    Zkv=tanh(f*Vkv);%función de activación nodos de salida
    Zv=round(Zkv);    
    ER=Zv-PVAL(3,n);
    h=h+1*(ER==0);
  
    clear y;clear z;clear zu;clear zv;clear ER;
    clear Netao
    clear Netas
end
h %cantidad de aciertos
