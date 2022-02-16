%MLP: XOR
%Juan Camilo Amaya
%Manuela Cerón Viveros
%Juan Daniel Muñoz

clc;
clear all;
close all;

%Entrada
X=[0 0 1 16;0 1 0 1];

%Target para XOR
tpk=[-1 1 1 -1];

%Inicio de los pesos
Wkj=zeros(1,2); 
Wji=zeros(2,2); 

%Defino capa oculta y de salida
Uj= zeros(2,1);
Yj= zeros(2,1);
Vk=0;
%Bias
bj=zeros(2,1);
bk=0;

%Factor de aprendizaje
E=0.01;

l=1;
o=0;
d=0;
e=0;

while o<4
    
    if l>4
        l=1;
    end

   %Capa oculta      
   Uj=Wji*X(:,l) + bj;
   Yj=tanh(Uj);
   
   %Capa de salida
   Vk=Wkj*Yj + bk;
   Zk=sign(purelin(Vk));
   Zk=Zk+(1*(Zk==0));
   
   %Error
   d=tpk(l)-Zk;
    
   diff1=1;
   diff2=(1-(tanh(Uj)).^2);
    
   o=(o+1)*(d==0);       
    
   %Actualización pesos capa de salida
   Wkj= Wkj+ (E*d.*diff1*Yj');
   bk= bk+(E*d.*diff1);
   
   %Actualización pesos capa oculta
    Wji=Wji+ (E*d.*diff1*Wkj*diff2*X(:,l))*ones(1,2);
    Wji=Wji+ (E*d.*diff1*Wkj*diff2*X(:,l))*ones(1,2);
    bj=bj + E*d.*diff2*diff1;
    
    l=l+1;
    e=e+1;
end


%%
%Validación

%Entrada
x1=1;
x2=1;

%Cálculo de la capa oculta
 Ujv=Wji(:,1)*x1 + Wji(:,2)*x2 + bj;
 
 Yjv=tanh(Ujv);
  
%Cálculo de la capa de salida 
Vkv=Wkj*Yjv + bk;
Zkv=sign(tanh(Vkv))

