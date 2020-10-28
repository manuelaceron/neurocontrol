%Simple Perceptron: Logic AND Function
%Manuela Cerón Viveros

clc;
clear all;
close all;

%Bias - Input 1 - Input 2
x=[1 1 1 1;0 0 1 1;0 1 0 1];

%Target
t=[-1 -1 -1 1];

N = size(t);
N = N(1);

%Weight initialization
w0=1.5;
w1=0.4;
w2=0.5;
W=[w0;w1;w2];

%Variable to iterate over x and t
k=1;
o=0;
delta=0;
e=0;
while o<N
    
    if k>N
        k=1;
    end
    
    %Signum function
    y=sign(W(1)*x(1,k)+W(2)*x(2,k)+W(3)*x(3,k));
    % y is 1 or -1
    y=y+(1*(y==0));
    % Delta calculation
    delta=t(k)-y;
    
    %Move to the next value
    o=(o+1)*(delta==0);
    
    %Weight calculation
    W=W+(1)*delta*(x(:,k));
    k=k+1;
    
    %Iterations counter
    e=e+1; 
end


%%
%Validation
%Ingress the input
x1=0;
x2=1;

%With the given weights, the model calculates the result for AND function
%The result will be 1 only when x1 and x2 are 1
yv=sign(W(1)+W(2)*x1+W(3)*x2)
