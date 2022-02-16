 %Octantes con salida lineal
clear all
clc;
close all;
display('Punto 1 - Tarea 2');
%Cantidad de iteraciones
N=1500;
Ne=input('Número de neuronas en la capa oculta ~4=');

X=zeros(3,N);
Tk=zeros(3,N);

for i=1:N
    %Generación de la entrada
    for j = 1:3
       %Generación de puntos en el rango [-10 -7] [-5 -2] [2 5] [7 10] 
       a=[3*rand(1,1)-10, 3*rand(1,1)-5, 3*rand(1,1)+2, 3*rand(1,1)+7];
       X(j,i)=a(randi(4));
    end
    
    %Target [-10 -7]=-2 [-5 -2]= -1 [2 5]=1 [7 10]=2
    for k= 1:3
      Tk(k,i)=(-2*(X(k,i)>=-10)*(X(k,i)<=-7))+(-1*(X(k,i)>=-5)*(X(k,i)<=-2))...
+(1*(X(k,i)>=2)*(X(k,i)<=5)) + (2*(X(k,i)>=7)*(X(k,i)<=10));
    end
end

%Inicialización de los pesos
Wji=-0.5+(0.5).*rand(Ne,3);%Entrada - capa oculta
Wkj=-0.5+(0.5).*rand(3,Ne);%Oculta - capa salida 


Tkx=Tk;
Xx=X;

for n=1:100;
    %Factor de aprendizaje variable
    alphax=1-(0.8*n)/100;
for k=1:N;
    %Cálculo de la capa oculta
    Uj=Wji*X(:,k);
    Yj=tanh(0.02*Uj); %Salida no lineal
    Dy=0.02*(1-Yj.^2);
    %Cálculo de capa de salida
    Vk=Wkj*Yj;
    Zk=purelin(Vk); %Salida lineal
    Dz=1;
    %Cálculo del error
    e(:,k)=Tk(:,k)-Zk;
    
    %cálculo del delta específico de salida
    deltaS=e(:,k).*Dz;
    
    %cálculo del delta específico oculto
    for j=1:Ne;
        deltaO(j)=Dy(j).*deltaS'*Wkj(:,j);
    end
    clear j
    
    %Actualización de los pesos oculta-salida
    for j=1:3;
        for i=1:Ne;
            Wkj(j,i)=Wkj(j,i)+alphax*deltaS(j).*Yj(i);
        end
        clear i
    end
    clear j
    
    %Actualización de los pesos entrada-oculta
    for j=1:Ne;
        for i=1:3;
            Wji(j,i)=Wji(j,i)+alphax*deltaO(j)*X(i,k);
        end
        clear i
    end
    clear j
    %Cálculo de error
    er(n,k)=0.5*sum(e(:,k).^2);
    clear Yj
    clear Zk
end
clear e
clear Uj
clear Vk
Xx=X(:,randperm(N));
Tkx=Tk(:,randperm(N));
end
%Calculo de error cuadrático medio
erx=(sum(er,2))/N;
figure(1)
plot(erx)

%%
%Validación manual

xp=input('coordenada en x=');
yp=input('coordenada en y=');
zp=input('coordenada en z=');
ve=[xp;yp;zp];
Netaop=Wji*ve;
Yj=tanh(0.02*Netaop);
Netasp=Wkj*Yj;
%Capa de salida
Zk=round(purelin(Netasp))



%%
%Validación completa

for i=1:100
    for j = 1:3
       av=[3*rand(1,1)-10, 3*rand(1,1)-5, 3*rand(1,1)+2, 3*rand(1,1)+7];
        Xval(j,i)=av(randi(4));
    end
   
    
    for k= 1:3
%      
    Tkv(k,i)=(-2*(Xval(k,i)>=-10)*(Xval(k,i)<=-7))+(-1*(Xval(k,i)>=-5)*(Xval(k,i)<=-2))...
    +(1*(Xval(k,i)>=2)*(Xval(k,i)<=5)) + (2*(Xval(k,i)>=7)*(Xval(k,i)<=10));
    end
end

h=0;
clear y;clear z;
for n=1:100;
    Ujv=Wji*Xval(:,n);%entrada neta en los nodos ocultos
    Yjv=tanh(0.02*Ujv);%función de transferencia nodos ocultos
    Vkv=Wkj*Yjv;%entrada neta en los nodos de salida
    Zkv=round(purelin(Vkv));
    Zv=Zkv;    

    ER=Zv-Tkv(:,n);
    h=h+1*(ER==0);
  

    clear y;clear z;clear zu;clear zv;clear ER;
    clear Netao
    clear Netas
end
display (h, 'aciertos de 100'); 



