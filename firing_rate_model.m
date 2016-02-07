% This was a firing rate model of interactions between these 3 cell types,
% as a first pass model.

clear

number_points = 100;
time_step = 0.001;      % 1 msec timestep
total_time = number_points*time_step;   
time = 0:time_step: total_time-time_step;
r_A = zeros(number_points,1);
r_B = zeros(number_points,1);
r_C = zeros(number_points,1);

i = poissrnd(100,number_points,1);


r_A0 = poissrnd(60,number_points,1);
r_B0 = poissrnd(5,number_points,1);
r_C0 = poissrnd(3.5,number_points,1);

r_A(1) = 60;
r_B(1) = 5;
r_C(1) = 3.5;

alpha_A = 2;
alpha_B = 2;
alpha_C = 3;
beta_B = 2;
beta_C = 2;


for n = 2:number_points
    r_A(n) = time_step*((-alpha_A*r_A(n-1)) + r_A0(n-1) + i(n-1) - beta_B*(r_B(n-1))) + r_A(n-1);
    r_B(n) = time_step*((-alpha_B*r_B(n-1)) + r_B0(n-1) + i(n-1)) + r_B(n-1);
    r_C(n) = time_step*((-alpha_C*r_C(n-1)) + r_C0(n-1) -beta_C*r_A(n-1)) + r_C(n-1);
end;

r_A = max(r_A, zeros(number_points,1));
r_B = max(r_B, zeros(number_points,1));
r_C = max(r_C, zeros(number_points,1));

figure(1)
clf
plot (time, r_A,'b')
hold on
plot(time,r_B,'m')
plot(time,r_C,'g')

A = fft(r_A);
B = fft(r_B);
C = fft(r_C);

figure(2)
subplot(3,1,1)
plot(A,'b')
subplot(3,1,2)
plot(B,'m')
subplot(3,1,3)
plot(C,'g')


