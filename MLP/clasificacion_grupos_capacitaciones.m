clc;
clear all;
close all;

% 4 grupos de capacitación, sin distinción de género:
% G1= empleado con estudios [0 0]
% G2= empleado sin estudios [0 1]
% G3= desempleado con estudios [1 0]
% G4= desempleado sin estudios [1 1]

X=[0 0 0 ; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1 ; 1 1 0; 1 1 1]';
tpk= [1 2 3 4 5 6 7 8];
ERROR=0;

n=2; %Número de neuronas en capa oculta
Wkj=zeros(1,n); 
Wji=zeros(n,3); 
bj=zeros(n,1);
bk=0;

Uj= zeros(n,1);
Yj= zeros(n,1);
Vk=0;
Zk=0;
E=0.1;

l=1;
o=0;
d=0;
e=0;

while o<8
    
    if l>8
        l=1;
    end
    
      
   %Capa oculta
   Uj=Wji*X(:,l) + bj; %Neta
   Yj=tanh(0.5*Uj); %Neta evaluada en función de activación
   dy=0.5*(1-(Yj).^2); 
   
   %Capa de salida
   Vk=Wkj*Yj + bk; %Neta
   Zk=round(purelin(Vk)); %Neta evaluada en función de activación
   dz=1;
   
   %Error
   d=tpk(l)-Zk;
   o=(o+1)*(d==0);
   
   %Delta capa de salida
   deltaS=d.*dz;
   
   %Delta capa oculta
   for i=1:n
    deltaO(i)=dy(i).*deltaS*Wkj(:,i);
   end
   clear i
    
   %Actualización pesos capa salida
    for i=1:n;
        Wkj(i)=Wkj(i)+E*deltaS.*Yj(i);
    end
    clear i
    bk= bk+(E*d*dz); %Bias

    %Actualización pesos capa oculta
    for j=1:n
        for i=1:3;
            Wji(j,i)=Wji(j,i)+E*deltaO(j)*X(i,l);
        end
        clear i
    end
    clear j
   bj=bj + E*d*dy*dz; %Bias
   
   l=l+1;
   e=e+1;
   
   
end


%%
%Encuestas
X=randi([0 1], 3, 10)
%Grupos de capacitación
g1=0;
g2=0;
g3=0;
g4=0;
g5=0;
g6=0;
g7=0;
g8=0;


for i= 1: length(X)
    
   %Uj=Wji*X(2:3,i) + bj;
   Uj=Wji*X(:,i) + bj;
   Yj=tanh(Uj);
   Vk=Wkj*Yj + bk;
   Zk=round(purelin(Vk));

   switch Zk
       case 1
       g1=g1+1;
       case 2
       g2=g2+1;
       case 3
           g3=g3+1;
       case 4
           g4=g4+1;
       case 5
           g5=g5+1;
       case 6
           g6=g6+1;
       case 7
           g7=g7+1;
       case 8
           g8=g8+1;
   end
end

res={'Grupo1', 'Grupo2', 'Grupo3', 'Grupo4', 'Grupo5', 'Grupo6', 'Grupo7', 'Grupo8';g1 g2 g3 g4 g5 g6 g7 g8}



