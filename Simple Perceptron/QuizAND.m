%Perceptrón simple: AND
%Juan Camilo Amaya
%Manuela Cerón Viveros
%Juan Daniel Muñoz

clc;
clear all;
close all;
%Bias - Entrada 1 - Entrada 2
x=[1 1 1 1;0 0 1 1;0 1 0 1];
%Target
t=[-1 -1 -1 1];
%Inicialización de los pesos
w0=1.5;
w1=0.4;
w2=0.5;
W=[w0;w1;w2];

k=1;
o=0;
d=0;
e=0;
while o<4
    
    if k>4
        k=1;
    end
    
    y=sign(W(1)*x(1,k)+W(2)*x(2,k)+W(3)*x(3,k));
    y=y+(1*(y==0))
    d=t(k)-y
    o=(o+1)*(d==0)
    W=W+(1)*d*(x(:,k))
    k=k+1
    
    e=e+1 %Contador de iteraciones
end


%%
%Validación
%Entrada
x1=1;
x2=1;

%Cálculo de la capa oculta
yv=sign(W(1)+W(2)*x1+W(3)*x2)
