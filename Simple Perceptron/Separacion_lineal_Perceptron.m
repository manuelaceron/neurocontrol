%Separación lineal: Perceptrón simple
%Juan Camilo Amaya
%Manuela Cerón Viveros
%Juan Daniel Muñoz

clear all;
close all;
clc
N=1000; %Número de muestras
r1=10; %Radio mayor
r2=6; %Radio menor
PX1=ones(4,N); %Media Luna Positiva
PX2=-1*ones(4,N); %MEdia luna negativa
hold on
k1=0; 
k2=0;
d=1; %Distancia de separación entre medias lunas
r=8; %SEparación horizontal entre medias lunas


%Generación de la entrada
while k1 <= N-1
    %Coordenadas entre x y y MEdia luna Positiva
    px1=randi([-100*r1,100*r1])/100;
    py1=randi([0,100*r1])/100;
    
    a=sqrt(px1^2+py1^2); %Magnitud
    
    %Rango entre r1 y r2
    if a<r1 && a>r2
        PX1(1,k1+1)=px1;
        PX1(2,k1+1)=py1;
        k1=k1+1
    end
end

while k2 <= N-1
    %Coordenadas entre x y y Media luna Negativa
    px2=randi([-100*r1,100*r1])/100;
    py2=(randi([-100*r1,0])/100);
    
    b=sqrt(px2^2+py2^2); %Magnitud
    
    %Rango entre r1 y r2
    if b<r1 && b>r2
        PX2(1,k2+1)=px2+r;
        PX2(2,k2+1)=py2-d;
        k2=k2+1;
    end
xlim([-20 20])
ylim([-20 20])
end


PX1(3,:)=PX1(3,:).*(-1); %Bias
P=[PX1 PX2]; %Media lunas positiva y negativa

figure(1)
hold on
plot(PX1(1,:),PX1(2,:),'.r')
hold on
plot(PX2(1,:),PX2(2,:),'.b')

%Inicialización de los pesos
w0=0.0;
w1=0.0;
w2=0.0;
W=[w0;w1;w2];
EPOCH=50;
NE=2000;

o=0;
d=0;
e=0;

%%
%Entrenamiento
for n=1:EPOCH
    for k=1:NE
    y=sign(W(1)*P(1,k)+W(2)*P(2,k)+W(3)*P(3,k));
    y=y+(1*(y==0));
    d=P(4,k)-y;
    o=(o+1)*(d==0);
    W=W+(1)*d*P(1:3,k);
    k=k+1;
    e=e+1;
    ERROR(n,k)=0.5*sum(d^2);
    end
    clear d
end
ERRORX=(sum(ERROR,2))/NE;
figure (2)
ylim([0 0.1])
plot(ERRORX); hold on;
%%
%Validación
Nv=2*N;
PX1v=ones(2,Nv);
PX2v=ones(2,Nv);
dv=1;
k3=0;
k4=0;
%Generación de entrada de Validación
while k3 <= Nv-1
    px1v=randi([-100*r1,100*r1])/100;
    py1v=randi([0,100*r1])/100;
    
    av=sqrt(px1v^2+py1v^2);
    
    if av<r1 && av>r2
        PX1v(1,k3+1)=px1v;
        PX1v(2,k3+1)=py1v;
        k3=k3+1
    end
end

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

Pv=[PX1v PX2v];

%Generación de la ecuación de la recta
for k=1:4000
m=-(W(1)/W(2));
b=(W(3)/W(2));
y1(k)=m* Pv(1,k)+b;
end

figure(3)
hold on
plot(PX1v(1,:),PX1v(2,:),'.r')
hold on
plot(PX2v(1,:),PX2v(2,:),'.b')
hold on
plot( Pv(1,:), y1)