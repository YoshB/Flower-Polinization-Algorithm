% Polinizacion de la flor


clc;
clear all;
close all;

% Choosing a function (Ej. 0, 1, 2... etc.)
f_objetivo = 2;

switch f_objetivo
    case 0
        %Griewank
       f = @(x,y) ((x.^2 + y.^2)/4000)-cos(x).*cos(y/sqrt(2)) + 1;
       U = [10 10];
       L = [-10 -10];
       %Minimum = 0; x= 0, y= 0
       
    case 1
        %Rastrigin
       f = @(x,y)  10*2 + x.^2 - 10 .*cos(pi*x)+ y.^2 - 10 .* cos(pi*y);
       U = [5 5];
       L = [-5 -5];
       %Minimum = 0; x = 0, y = 0 
       
      
    case 2
         %DropWave
    f = @(x,y) - ((1 + cos(12*sqrt(x.^2+y.^2))) ./ (0.5 *(x.^2+y.^2) + 2));
    U = [2 2];
    L = [-2 -2];
       %Minimum = -1; x = 0, y =0
      
    case 3
        %Esphere
        f = @(x,y) (x-2).^2 + (y-2).^2;
        U = [5 5];
        L = [-5 -5];
        %Minimum = 0; x = 2, y = 2 
       
    otherwise
        disp("Choose a valid value for the function")
        return
end

% % % % % % % % % % % % % % % % % % % % % % % % % % 
%contour of the function
[X,Y] = meshgrid(L(1):0.2:U(1), L(2):0.2:U(2));
Z = f(X,Y);
grid on;
contour(X,Y,Z,40);
hold on
% % % % % % % % % % % % % % % % % % % % % % % % % % 

%Number of variables of the problem (X,Y) in this case
D = 2;
%Population size
N = 50;
%Iterations
G = 150;

%Probability critery
p= 0.8;

% Calculating the value for the Lévy fly
lambda = 1.5;
sigma2= gamma(1+lambda) / (lambda * gamma((1+lambda)/2));
sigma2 = sigma2 * ( sin(pi*lambda/2) / (2^((lambda-1)/2)) );
simga2 = sigma2^(1/lambda);

%list with the population
individuos = zeros(D,N);


% Inicialization
r = rand(D,N);
for i=1:D
    individuos(i,:) = L(i) + (U(i) - L(i)) .* r(i,:);
end

fitness = f(individuos(1,:),individuos(2,:));

for i=1:G
    %Choose the better solution
    [~,value] = min(fitness);
    best = individuos(:,value);
    
    for j=1: N
        
        r = rand();
        
        if r < p
            %Calculating 'L' for the Lévy fly
            u = normrnd(0,sigma2,[D,1]);
            v = normrnd(0,1,[D,1]);
            
            Lev = u./ ( abs(v).^(1/lambda));
            
            y = individuos(:,j) + Lev.*(individuos(:,value)-individuos(:,j));
            
        else
            rn = rand();
            %Choosing three individuals
            k = j;
            n = j;
            while k == j
                k = randi([1,N]);          
            end
            while n == k || n == j
                n = randi([1,N]); 
            end
            
            y = individuos(:,j) + rn*(individuos(:,k) -individuos(:,n));
        end
        
        if f(y(1),y(2)) < fitness(j)
            individuos(:,j) = y;
        end
    end
    
    fitness = f(individuos(1,:),individuos(2,:));

    
%     Graphics
 axis([L(1) U(1) L(2) U(2)]);
 h = plot(individuos(1,:),individuos(2,:), 'rx');
 pause(0.1);
 delete(h);
    
end

[fit_best, value] = min(fitness);
best = individuos(:,value)
fit_best


