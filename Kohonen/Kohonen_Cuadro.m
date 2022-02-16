%Códido por David Ramírez
%Adaptado por Juan Camilo Amaya
%Manuela Cerón Viveros
%Juan Daniel Muñoz
%RED KOHONEN

%SOM_número de neuronas igual que número de puntos
%SOM_posición de neuronas predeterminada_vecindad por matriz de posición
close all;
clc;
display('Kohonen');
%*************** Generación de datos en el especio sensorial/entrada
nump=input('Número de datos: nump= ');%número total de datos
neux=input('Númerode posiciones x en la capa neural: neux=');
neuy=input('Númerode posiciones y en la capa neural: neuy=');
entrfx=zeros(1,nump);
entrfy=zeros(1,nump);
entrfz=zeros(1,nump);
entrfu=zeros(1,nump);
%***************
    for i=1:nump;
        test=rand(1);
        entrfx(1,i)=(test<=0.5)*(-10+5*rand(1))+(test>0.5)*(5+5*rand(1));
        test=rand(1);
        entrfy(1,i)=(test<=0.5)*(-8+5*rand(1))+(test>0.5)*(3+5*rand(1));
        test=rand(1);
        entrfz(1,i)=(test<=0.5)*(-7+1*rand(1))+(test>0.5)*(6+1*rand(1));
        test=rand(1);
        entrfu(1,i)=(test<=0.5)*(-10+2*rand(1))+(test>0.5)*(8+2*rand(1));
    end
    clear i
%*******************
redpun=[entrfx;entrfy;entrfz;entrfu];%matriz de datos
%****************
indexs=randperm(nump);
redpunp=redpun(:,indexs);%permutación en redneu
redpunpt=redpunp;
for i=1:99;%concatenación
    redpunpt=[redpunpt redpun(:,randperm(nump))];
end
clear i
%***************
neutot=neux*neuy;%número total de neuronas
entrdx=zeros(1,neutot);
entrdy=zeros(1,neutot);
for i=1:neutot;
    entrdx(1,i)=ceil(i/neuy);
    entrdy(1,i)=mod(i,neuy)*(mod(i,neuy)~=0)+(neuy)*(mod(i,neuy)==0); 
end
clear i
redneu=[entrdx;entrdy];%posiciones de las neuronas
%**************** Inicialización de pesos
pesosx=-0.5+rand(1,neutot);
pesosy=-0.5+rand(1,neutot);
pesosz=-0.5+rand(1,neutot);
pesosu=-0.5+rand(1,neutot);
pesosi=[pesosx;pesosy;pesosz;pesosu];%pesos sinápticos iniciales
%****************
Rxi=ceil(1*min(neux,neuy));
Ryi=ceil(1*min(neux,neuy));
alphai=0.9;
alphaf=0.1;
Rf=1;
%*******************
pesos=pesosi;
QQ=0;
for i=1:length(redpunpt);
    for j=1:neutot;
        normaresta(j)=norm(redpunpt(:,i)-pesos(:,j));
    end
    clear j
    Rx(i)=Rxi+(Rf-Rxi)*((i-1)/length(redpunpt));
    Ry(i)=Ryi+(Rf-Ryi)*((i-1)/length(redpunpt));
    alpha(i)=alphai+(alphaf-alphai)*((i-1)/length(redpunpt));
    winner(i)=find(normaresta==min(normaresta),1);
    pesos(:,winner(i))=pesos(:,winner(i))+ alpha(i)*(redpunpt(:,i)-pesos(:,winner(i)));
    winnerx(i)=ceil(winner(i)/neuy);
    winnery(i)=mod(winner(i),neuy)*(mod(winner(i),neuy)~=0)+(neuy)*(mod(winner(i),neuy)==0);
    for j=1:neutot;
        if abs(ceil(j/neuy)-winnerx(i))<=Rx(i)&&...
                abs((mod(j,neuy)*(mod(j,neuy)~=0)+(neuy)*(mod(j,neuy)==0))-winnery(i))<=Ry(i);
            pesos(:,j)=pesos(:,j)+ alpha(i)*(redpunpt(:,i)-pesos(:,j));
        end
    end
    QQ=QQ+1;
    if QQ==round((length(redpunpt))/3);
        display('1 tercio');
    elseif QQ==round(2*(length(redpunpt))/3);
        display('2 tercios');
    elseif QQ==length(redpunpt);
        display('final');
    end
end
clear i
%****************
%Definir rangos para colores por categoria
funct=@(x)(99/15)*x+(1-99/15);
for i=1:neutot
    if pesos(1,i)>0 && pesos(2,i)>0 && pesos(3,i)>0 && pesos(4,i)>0;
        map(redneu(1,i),redneu(2,i))=funct(16);
    elseif pesos(1,i)>0 && pesos(2,i)>0 && pesos(3,i)>0 && pesos(4,i)<0;
        map(redneu(1,i),redneu(2,i))=funct(15);
    elseif pesos(1,i)>0 && pesos(2,i)>0 && pesos(3,i)<0 && pesos(4,i)>0;
        map(redneu(1,i),redneu(2,i))=funct(14);
    elseif pesos(1,i)>0 && pesos(2,i)<0 && pesos(3,i)>0 && pesos(4,i)>0;
        map(redneu(1,i),redneu(2,i))=funct(13);
    elseif pesos(1,i)<0 && pesos(2,i)>0 && pesos(3,i)>0 && pesos(4,i)>0;
        map(redneu(1,i),redneu(2,i))=funct(12);
    elseif pesos(1,i)>0 && pesos(2,i)>0 && pesos(3,i)<0 && pesos(4,i)<0;
        map(redneu(1,i),redneu(2,i))=funct(11);
    elseif pesos(1,i)>0 && pesos(2,i)<0 && pesos(3,i)<0 && pesos(4,i)>0;
        map(redneu(1,i),redneu(2,i))=funct(10);
    elseif pesos(1,i)<0 && pesos(2,i)<0 && pesos(3,i)>0 && pesos(4,i)>0;
        map(redneu(1,i),redneu(2,i))=funct(9);
    elseif pesos(1,i)>0 && pesos(2,i)<0 && pesos(3,i)>0 && pesos(4,i)<0;
        map(redneu(1,i),redneu(2,i))=funct(8);
    elseif pesos(1,i)<0 && pesos(2,i)>0 && pesos(3,i)<0 && pesos(4,i)>0;
        map(redneu(1,i),redneu(2,i))=funct(7);
    elseif pesos(1,i)<0 && pesos(2,i)>0 && pesos(3,i)>0 && pesos(4,i)<0;
        map(redneu(1,i),redneu(2,i))=funct(6);
    elseif pesos(1,i)>0 && pesos(2,i)<0 && pesos(3,i)<0 && pesos(4,i)<0;
        map(redneu(1,i),redneu(2,i))=funct(5);
    elseif pesos(1,i)<0 && pesos(2,i)<0 && pesos(3,i)<0 && pesos(4,i)>0;
        map(redneu(1,i),redneu(2,i))=funct(4);
    elseif pesos(1,i)<0 && pesos(2,i)>0 && pesos(3,i)<0 && pesos(4,i)<0;
        map(redneu(1,i),redneu(2,i))=funct(3);
    elseif pesos(1,i)<0 && pesos(2,i)<0 && pesos(3,i)>0 && pesos(4,i)<0;
        map(redneu(1,i),redneu(2,i))=funct(2);
    elseif pesos(1,i)<0 && pesos(2,i)<0 && pesos(3,i)<0 && pesos(4,i)<0;
        map(redneu(1,i),redneu(2,i))=funct(1);
    end
end
clear i
%*************************
figure(1)      
surf(map)
        
 

    
    
    
    
        
        
        