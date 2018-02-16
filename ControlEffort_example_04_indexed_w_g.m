% This file has been added to the git repository for test-only purpose.
% This example maps the stabilizing set indexed by a gain crossover 
% frequency to a 2-plane in which x-axis is the
% H-infinity norm of the error transfer function and y-axis is the
% H-infinity norm of the control transfer function.

% initialize
close all;
clear;
tic
% plant
syms s;
N = [1 -2];         % Numerador
D = [1 4 3];        % Denominador
P_tf=tf(N,D);
N_s=poly2sym(N,s);
D_s=poly2sym(D,s);
P=N_s/D_s;
theta_rad=-pi:0.01:pi;

% load stabilizing set
load stab_pi_plant_00.mat;
[xv, yv] = poly2cw([Ki_bounds(:,1)' flipud(Ki_bounds(:,2))'],[Kp_range' flipud(Kp_range)']);

w_g_vector = linspace(0.1, 3, 30);
test_points_cell = cell(numel(w_g_vector),1);
for idx=1:numel(w_g_vector)
    w_g = w_g_vector(idx);
    Pr  = real(subs(P,s,1i*w_g));
    Pi  = imag(subs(P,s,1i*w_g));
    bb  = 1/(Pr^2 + Pi^2);
    aa  = (w_g^2)*bb;
    xq_g = double(sqrt(aa)*cos(theta_rad));
    yq_g = double(sqrt(bb)*sin(theta_rad));
    [in, on] = inpolygon(xq_g, yq_g, xv, yv);
    test_points_cell(idx)={[xq_g(in&~on);yq_g(in&~on)]};
end

mapping_points_cell = cell(numel(w_g_vector),1);
figure; hold on;
for idx = 1:numel(w_g_vector)
    test_points_matrix = cell2mat(test_points_cell(idx));
    x=test_points_matrix(1,:)';
    y=test_points_matrix(2,:)';
    reru_matrix = NaN(2,numel(x));
    for idy = 1:numel(x)
        Ki = x(idy);
        Kp = y(idy);
        C = pid(Kp,Ki);
        sys_error = minreal(1/(1+P_tf*C));
        sys_control = minreal(C/(1+P_tf*C));
        [re, f_re] = getPeakGain(sys_error);
        [ru, f_ru] = getPeakGain(sys_control);
        reru_matrix(:,idy) = [re; ru];
    end
    plot(reru_matrix(1,:),reru_matrix(2,:),'k-x');
end
hold off;
% figure;
% hold on;
% fill([Ki_bounds(:,1)' flipud(Ki_bounds(:,2))'],[Kp_range' flipud(Kp_range)'],'y');
% %plot(xq_g,yq_g,'bx');
% %plot(xq_g(in&~on),yq_g(in&~on),'rx');
% xlabel('K_I');
% ylabel('K_P');
% axis([-4.4 1.2 -4.5 2.0]);
% hold off;
% toc
