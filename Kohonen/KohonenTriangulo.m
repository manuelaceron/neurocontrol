%Códido por David Ramírez
%Entrada triangular adaptada por:
%Juan Camilo Amaya
%Manuela Cerón Viveros
%Juan Daniel Muñoz
%RED KOHONEN


%SOM_posición de neuronas predeterminada_vecindad por cardinal
clear all;close all;clc;
display('Kohonen Triángulo');
numpy=input('Nro de posiciones y en el espacio físico: numpy=');
numpx=input('Nro de posiciones x en el espacio físico: numpx=');
neuy=numpy/2; %Ubicar en el centro
neux=numpx/2;
neutot=input('Nro total de neuronas: neutot=');
nump=neutot;
Np=input('Nro de veces que se presentan las entradas Np=');
Rxi=0.5*ceil(1*min(neux,neuy)); %vecidad?
alphai=0.9;
alphaf=0.1;
Rf=0;
%****
r=1;
k=0;
%Organización de las neuronas en circulo
while k <= neutot-1
    px=randi([-100*r,100*r])/100;
    py=randi([-100*r,100*r])/100;
    if sqrt(px^2+py^2)<r
         pesosi(1,k+1)=px+neuy;
         pesosi(2,k+1)=py+neuy;
        k=k+1;
    end
end
%*****
entrdx=zeros(1,nump);
entrdy=zeros(1,nump);

%Organización del triangulo

x1=1:(neuy-1)/numpy:neuy-(neuy-1)/numpy;
x2=neuy:(neuy-1)/numpy:numpy-1-(neuy-1)/numpy;
x3=1:(numpy-2)/numpy:(numpy-1)-(numpy-2)/numpy;

 for j=1:numpy
y1(j)=2*x1(j);
y2(j)=-2*x2(j)+2*numpy;
 end
y3= ones(1,numpy)+1;
redpun=[x1 x2 x3; y1 y2 y3];
 
figure(1)
plot(redpun(1,:),redpun(2,:),'r*'); 
hold on
plot(pesosi(1,:),pesosi(2,:),'b*')


redpunpt=redpun;
pesos=pesosi;
for i=1:Np*neutot;
    indexs=randperm(nump); %Revolver
    redpunpt=redpun(:,indexs);
    
     k=mod(i,neutot)*(mod(i,neutot)~=0)+(neutot)*(mod(i,neutot)==0);
    for j=1:neutot;
        normaresta(j)=norm(redpunpt(:,k)-pesos(:,j));
    end
    clear j
    Rx(i)=Rxi+(Rf-Rxi)*((i-1)/(Np*neutot));
    alpha(i)=alphai+(alphaf-alphai)*((i-1)/(Np*neutot));
    winner(i)=find(normaresta==min(normaresta),1);
    pesos(:,winner(i))=pesos(:,winner(i))+ alpha(i)*(redpunpt(:,k)-pesos(:,winner(i)));
    for j=1:neutot;
        if abs(winner(i)-j)<=Rx(i);
           pesos(:,j)=pesos(:,j)+ alpha(i)*(redpunpt(:,k)-pesos(:,j));
%            Graficación del cambio en los pesos en tiempo real
%            anim=zeros(2,neutot);
%            anim(:,j)=pesos(:,j); 
%            plot(anim(1,j),anim(2,j),'b*'); 
%            M(j)=getframe;
        end
    end
 i
end
% movie(M)
figure(2)
plot(redpun(1,:),redpun(2,:),'r*'); 
hold on
plot(pesos(1,:),pesos(2,:),'b*')
hold off
figure(3)
plot(pesos(1,:),pesos(2,:),'b*');    