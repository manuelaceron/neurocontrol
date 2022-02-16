%Juan Camilo Amaya
%Juan Daniel Muñoz
%Manuela Cerón 
%Red Neuronal RBF

clear all;
close all;
clc;

display('Entrenamiento de centroides RBF');
np=input('Cantidad de Neuronas:'); %Cantidad de patrones totales
n4=round(np/4); %patrones para cada círculo 
nc=input('Cantidad de Centroides:');%cantidad de centroides 
error=0;% error admitido para el entrenamiento de los centroides 


%Vector de entrada: datos aleatorios con 4 centroides dentro de un radio de 2
%****** Generación de patrones para cada círculo ******
k=0;
k1=0;
k2=0;
k3=0;
r=2;

P1=ones(2,n4);
P2=ones(2,n4);
P3=ones(2,n4);
P4=ones(2,n4);

while k <= n4-1
    px=(randi([-100*r,100*r])/100);
    py=(randi([-100*r,100*r])/100);
    px_1=px+2;
    py_1=py+2;
    if sqrt(px^2+py^2)<r
        P1(1,k+1)=px_1;
        P1(2,k+1)=py_1;
        k=k+1;
       
    end
end   

while k1 <= n4-1
    
    px1=(randi([-100*r,100*r])/100);
    py1=(randi([-100*r,100*r])/100);
    px_2=px1-2;
    py_2=py1+2;
        if sqrt(px1^2+py1^2)<r
        
        P2(1,k1+1)=px_2;
        P2(2,k1+1)=py_2;
        k1=k1+1;
        end
end

while k2 <= n4-1
    
    px2=(randi([-100*r,100*r])/100);
    py2=(randi([-100*r,100*r])/100);
    px_3=px2+2;
    py_3=py2-2;
    if sqrt(px2^2+py2^2)<r
        
        P3(1,k2+1)=px_3;
        P3(2,k2+1)=py_3;
        k2=k2+1;
    end
end


while k3 <= n4-1
    
    px3=(randi([-100*r,100*r])/100);
    py3=(randi([-100*r,100*r])/100);
    px_4=px3-2;
    py_4=py3-2;
    if sqrt(px3^2+py3^2)<r
       
        P4(1,k3+1)=px_4;
        P4(2,k3+1)=py_4;
        k3=k3+1;
    end
end
%******************************************************************

P=[P1(1,:) P2(1,:) P3(1,:) P4(1,:); P1(2,:) P2(2,:) P3(2,:) P4(2,:)]; %vector total de patrones
figure(1)
plot(P(1,:),P(2,:),'.b')
hold on;
%Declaración de centroides iniciales

ci=P(:,randperm(np,nc));% elección aleatorea de centroides iniciales
plot(ci(1,:),ci(2,:),'*k','linewidth',1.5) %Centroides iniciales

dis=zeros(np,nc);
%%

cen=zeros(np,1);
o=0;
iteraciones=0;
promxy_anterior=zeros(2,nc);
total=zeros(1,nc);
promxy=ci;

while o<nc %ciclo que termina cuando los centroides centroide(i) = centroide(i) (que no actualizen mas)
winx=zeros(np,nc);
winy=zeros(np,nc);
for i=1:np %ciclo de 1 a la cantidad de patrones totales
    
    for j=1:nc %ciclo de 1 a cantidad de centroides
        dis(i,j)=sqrt(sum((P(:,i)-promxy(:,j)).^2)); %comparación de posición de los patrones con cada centroide
        pos=find(dis(i,:)==min(dis(i,:)),1); %centroide con el patron más cercano 
         
        cen(i,1)=pos; %se guardan las posiciones en un vector 
        
    end
        
        %*** guarda cada valor X y Y de los patrones en columnas que representan a cada centroide 
        winx(i,pos)=P(1,i); 
        winy(i,pos)=P(2,i);
        %***
        
end

clear i;

for i=1:nc %cilco de 1 a cantidad de centroides
        
        
        %*** encuentra la cantidad de patrones pertenecientes a cada patron
        Re=find(cen==i); 
        total(i)=length(Re);
        %***
        
        %*** promedio de los centroides cercanos a cada patron 
        promxy(1,i)=sum(winx(:,i))/total(i);
        promxy(2,i)=sum(winy(:,i))/total(i);
        %***
        
        %*** comparación de valores de centroides del ciclo anterior con el actual
        o=o+1*(abs(promxy_anterior(1,i)-promxy(1,i))<=error)*(abs(promxy_anterior(2,i)-promxy(2,i))<=error);
        o=o*(abs(promxy_anterior(1,i)-promxy(1,i))<=error)*(abs(promxy_anterior(2,i)-promxy(2,i))<=error);
        %***
end
        
        promxy_anterior=promxy; %actualización del promedio del ciclo anterior
        iteraciones=iteraciones+1; %conteo de iteraciones 
end
 plot(promxy(1,:),promxy(2,:),'ok','linewidth',1.5); grid on;

clear i j;

%*** Calculo de sigma por el método 1
for i=1:nc
    for j=1:nc
       dis_cent(i,j)=sqrt(sum((promxy(:,i)-promxy(:,j)).^2));
    end
end
 sigma1=sum(dis_cent)/nc;
 %***
 
 %*** calculo de sigma por el método 2
 sigma2=sum(dis)/np;
 %***
 
 %% Capa de salida 
 % la capa de salida para la RBF se compone de un método de caracterización
 % entrenada donde se conoce las categorias a clasificar, por tal razón
 % éste método no hace uso de pruning.
 
 %Para el uso de las mismas variables ya abarcadas en el trenamiento
 %anterior se limpian todas
clear i j k k1 k2 k3 P1 P2 P3 P4 y x z ;

%MLP sin momentum/Conjunto de entrenamiento con repetición
display('Capa de salida RBF');
NE=input('Número de vectores de entrada de entrenamiento (~2000): NE='); % cantidad de vectores de entrenamiento
Nodos=input('Número de neurnas de salida: Nodos='); %nodos capa de salida

%se usa la forma de la entrada del entrenamiento de centroides ya que es el
%misma entrada para todo el sistema 

%************************
n4=round(NE/4);

k=0;
k1=0;
k2=0;
k3=0;
r=2;

P1=ones(2,n4);
P2=ones(2,n4);
P3=ones(2,n4);
P4=ones(2,n4);

while k <= n4-1
    px=(randi([-100*r,100*r])/100);
    py=(randi([-100*r,100*r])/100);
    px_1=px+2;
    py_1=py+2;
    if sqrt(px^2+py^2)<r
        P1(1,k+1)=px_1;
        P1(2,k+1)=py_1;
        k=k+1;
       
    end
end   

while k1 <= n4-1
    
    px1=(randi([-100*r,100*r])/100);
    py1=(randi([-100*r,100*r])/100);
    px_2=px1-2;
    py_2=py1+2;
        if sqrt(px1^2+py1^2)<r
        
        P2(1,k1+1)=px_2;
        P2(2,k1+1)=py_2;
        k1=k1+1;
        end
end

while k2 <= n4-1
    
    px2=(randi([-100*r,100*r])/100);
    py2=(randi([-100*r,100*r])/100);
    px_3=px2+2;
    py_3=py2-2;
    if sqrt(px2^2+py2^2)<r
        
        P3(1,k2+1)=px_3;
        P3(2,k2+1)=py_3;
        k2=k2+1;
    end
end


while k3 <= n4-1
    
    px3=(randi([-100*r,100*r])/100);
    py3=(randi([-100*r,100*r])/100);
    px_4=px3-2;
    py_4=py3-2;
    if sqrt(px3^2+py3^2)<r
       
        P4(1,k3+1)=px_4;
        P4(2,k3+1)=py_4;
        k3=k3+1;
    end
end
%************************

P=zeros(nc,NE);
P=[P1(1,:) P2(1,:) P3(1,:) P4(1,:); P1(2,:) P2(2,:) P3(2,:) P4(2,:)];
figure(2)
plot(P(1,:),P(2,:),'.b')
hold on;

entrenamiento=(-1)*ones(NE,nc);
dis2=zeros(NE,nc);
cen2=zeros(NE,1);

for i=1:NE
    for j=1:nc
        %se calcula la distancias de la entrada con cada centroide y se compara con el umbral de clasificación de cercanía 
        dis2(i,j)=(sqrt(sum((P(:,i)-promxy(:,j)).^2)))<0.8;
        %La salida de dis2 es binaria por lo que hay que encontar la posición del mayor 
        pos2=find(dis2(i,:)==max(dis2(i,:)),1)*((max(dis2(i,:)))>0);
        cen2(i,1)=pos2; % se guardan las posiciones en una variable
    
    end
    
    clear j;

end

 M=[dis2';cen2']; %vector de entrada y target 

%******************************
% Weoi=-0.25+(0.5).*rand(Nodos,nc);%matriz de pesos sinápticos entrada - capa oculta
% Wosi=-0.25+(0.5).*rand(1,Nodos);%matriz de pesos sinápticos capa oculta - capa salida 

Weoi=zeros(Nodos,nc);%matriz de pesos sinápticos entrada - capa oculta
Wosi=zeros(1,Nodos);%matriz de pesos sinápticos capa oculta - capa salida 
Wji=Weoi; Wkj=Wosi; 
%******************************
error=zeros(nc,NE);
Mx=M;
entr=cen2';
bj=zeros(Nodos,1);
bk=0;

Uj= zeros(Nodos,1);
Yj= zeros(Nodos,1);
Vk=0;
Zk=0;
E=0.1;
for n=1:100;
for k=1:NE;
   
   %Capa oculta
   Uj=Wji*Mx(1:nc,k) + bj; %Neta
   Yj=tanh(0.5*Uj); %Neta evaluada en función de activación
   dy=0.5*(1-(Yj).^2); 
   
   %Capa de salida
   Vk=Wkj*Yj + bk; %Neta
   Zk=round(purelin(Vk)); %Neta evaluada en función de activación
   %Zk=(tanh(Vk));
   dz=1;
   
   %Error
   d=Mx([nc+1],k)-Zk;
   o=(o+1)*(d==0);
   
   %Delta capa de salida
   deltaS=d.*dz;
   clear i
   %Delta capa oculta
   for i=1:Nodos
    deltaO(i)=dy(i).*deltaS*Wkj(:,i);
   end
   clear i
    
   %Actualización pesos capa salida
    for i=1:Nodos;
        Wkj(i)=Wkj(i)+E*deltaS.*Yj(i);
    end
    clear i
    bk= bk+E*d*dz; %Bias

    %Actualización pesos capa oculta
    
    for j=1:Nodos
        for i=1:nc;
            Wji(j,i)=Wji(j,i)+E*deltaO(j)*Mx(i,k);
        end
        clear i
    end
    
    clear j
    bj=bj + E*d*dy*dz; %Bias
    %**************************
    clear y
    clear z
   
end
clear error
clear Netao
clear Netas
Mx=M(:,randperm(NE));
%entr=entrenamiento(:,randperm(NE));
end
%*****************
%% Validación
%Se realiza una validación individual donde se busca clasificar la entrada
%ante la en centroide más cecano.

xp=input('xp=');
yp=input('yp=');
val=[xp yp]';
meta=(-1)*ones(1,nc);
dis3=zeros(1,nc);

for j=1:nc
    
        dis3(1,j)=(sqrt(sum((val(:,1)-promxy(:,j)).^2)))<0.8;
        pos3=(find(dis3(1,:)==max(dis3(1,:)),1))*((max(dis3(1,:)))>0);
       
        cen3(1,1)=pos3;

 end


ve=dis3';
Netaop=Wji(1,:)*dis3(1,:)' +bj;
y=tanh(0.5*Netaop);
Netasp=Wkj(1,1).*y(1,1) +bk;
z=round(purelin(Netasp)) %salida de la validación 









