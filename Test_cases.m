clear
clc
%%  Check that the 1D Minimization Algorithms by comparing the result with Rosenbrock s parabolic functon Graph
fprintf('-------------Applaying 1D minimization algorithms on Rosenbrock-----------\n');
f_Rose=Rosenbrock(-1.2,1.0);
a=0;
b=0.1;
t0=0.01;
error=0.01;
fib(a,b,error,f_Rose,'Min'); % The Fibonacci method
golden(a,b,error,f_Rose,'Min'); % The golden section method
Quad(t0,error,f_Rose); % The quadratic interpolation method
Cubic(t0,error,f_Rose); % The cubic interpolation method the most accurate methods of the above1
%%  Check that the 1D Minimization Algorithms by comparing the result with Powell functon Graph
fprintf('-------------Applaying 1D minimization algorithms on Powell problem------------\n');
f_powell=Powell(3,-1,0,1);
a=0;
b=0.1;
t0=0.01;
error=0.01;
fib(a,b,error,f_powell,'Min'); % The Fibonacci method
golden(a,b,error,f_powell,'Min'); % The golden section method
Quad(t0,error,f_powell); % The quadratic interpolation method
Cubic(t0,error,f_powell); % The cubic interpolation method the most accurate methods of the above1
%% Fletcher-Reeves CG Method
fprintf('-----------Applaying Fletcher-Reeves CG algorithms on Rosenbrock’s and Powell problem---------\n');
eR=0.1;
eP=0.1;
RX0=[-1.2;1];
PX0=[3;-1;0;1];
Fletcher_Reeves(eR,eP,RX0,PX0);
%% plotting compartion
lambda_range = linspace(-2, 2, 100);
% Evaluate the functions at different values of lambda
f_powell_values = f_powell(lambda_range);
f_Rose_values = f_Rose(lambda_range);
figure;
% Plot for Rosenbrock function
subplot(2, 1, 1);
plot(lambda_range, f_Rose_values, 'LineWidth', 2);
title('Rosenbrock Function along a Direction');
xlabel('\lambda');
ylabel('f(\lambda)');
grid on;
% Plot for Powell function
subplot(2, 1, 2);
plot(lambda_range, f_powell_values, 'LineWidth', 2);
title('Powell Function along a Direction');
xlabel('\lambda');
ylabel('f(\lambda)');
grid on;
%  the layout for better visualization
sgtitle('1D Minimization Algorithms on Powell and Rosenbrock Functions'); 
