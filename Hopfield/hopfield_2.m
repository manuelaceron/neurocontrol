
%Red Hopfield 
%Juan Camilo Amaya
%Manuela Cerón
%Juan Daniel Muñoz

clc; clear all;
display('Red Hopfield') 

%Cargo un data set de imagenes numéricas de matlab
 load('digit_train.mat', 'xTrainImages');

 %Patrones para aprendizaje cuadrados de igual dimensión
 uno=(xTrainImages{12});
 nueve=(xTrainImages{2});
 
 na=2; %cantidad de patrones para aprendizaje
 nv=4; %cantidad de patrones para validar
 
 %se convierten las imágenes a binario
 in1=im2bw(uno);
 in2=im2bw(nueve);
 
 in3=in1;
 in4=in2;
 
 %Dimensión
 [fm,cm]=size(in1); 
 
 %Distorsión de los patrones in3 e in4
  for i=16:20
    for j=10:15
      in3(i,j)=1;
      in4(i,j)=1;   
    end
  end 
  for i=3:10
    for j=13:20
      in3(i,j)=0;
      in4(i,j)=0;   
    end
  end 
 
 %Convierto matriz de patrones en un vector fila
 x(1,:)=reshape(in1, 1, []);
 x(2,:)=reshape(in2, 1, []);
 x(3,:)=reshape(in3, 1, []);
 x(4,:)=reshape(in4, 1, []);
 
 %Dimensión de la fila
 [fc,cc]=size(x);
 
 %Convierto patrones a -1 y 1
 for i=1:nv
    for j=1:cc
      if x(i,j)==0
          y(i,j)=-1;
          else
          y(i,j)=1;
      end
    end
 end  

%Aprendizaje de los dos primeros patrones 
for i =1:1:cc
    for j=1:1:cc
        p=0;
        if i ~= j
            for k=1:na 
                p=p+y(k,i).*y(k,j);
            end
        end
        pesos(i,j)=p;
    end
end

   clear i j k
   
   %Validación de los dos segundos patrones
    for k=na+1:nv
     for j=1:cc
        for i=1:cc
           if i ~= j
           h(i,j,k)=y(k,i).*pesos(j,i);
           end
       end
   end
   end
sumas=sum(h);
sig=sign(sumas);


 %Convierto a 0 y 1 para mostrar imagen binaria
 for i=1:nv
    for j=1:cc
      if sig(1,j,i)==-1
          out(i,j)=0;
          else
          out(i,j)=1;
      end
    end
 end
 
%Convierto vector fila a matriz cuadrada
val1=reshape(out(1,:), cm,[]);
val2=reshape(out(2,:), cm,[]);
val3=reshape(out(3,:), cm,[]);
val4=reshape(out(4,:), cm,[]);

%Imprimo
%Patrones aprendidos
subplot(2,4,1)
imshow(in1)
subplot(2,4,2)
imshow(in2)
%Patrones distorsionados para validación
subplot(2,4,3)
imshow(in3)
subplot(2,4,4)
imshow(in4)
%Salida de la red
subplot(2,4,6)
imshow(val3)
subplot(2,4,7)
imshow(val4)
