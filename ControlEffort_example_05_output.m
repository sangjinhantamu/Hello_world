% H-infinity norm trade off
% error function versus control function
% indexed by output function

close all;
clear;
syms s;

tic
% plant
N = [1 -2];         % Numerador
D = [1 4 3];        % Denominador
P_tf=tf(N,D);
N_s=poly2sym(N,s);
D_s=poly2sym(D,s);
P=N_s/D_s;
theta_rad=-pi:0.01:pi;

load stab_pi_plant_00.mat;
stabset_x = Ki_bounds(:,1)';
stabset_y = Kp_range';
[stabset_x, stabset_y] = poly2cw(double(stabset_x), double(stabset_y));

gamma_y = [1.1:0.1:1.9 2:16];
w_freqs=0.01:0.01:5;
test_points_cell = cell(numel(gamma_y),1);
for idx = 1:numel(gamma_y)
    gamma = gamma_y(idx);
    for idy = 1:numel(w_freqs)
        w=w_freqs(idy);
        Pr=real(subs(P,s,1i*w));
        Pi=imag(subs(P,s,1i*w));
        
        % ellipse for the output function PC/(1+PC)
        % (x - c1)^2 / a^2 + (y - c2)^2 / b^2 = 1
        c2 = (gamma^2*Pr)/((1 - gamma^2)*(Pr^2 + Pi^2));
        c1 = (gamma^2*w*Pi)/((1 - gamma^2)*(Pr^2 + Pi^2));
        bb = (gamma^2)/((1 - gamma^2)^2*(Pr^2 + Pi^2));
        aa = (w^2)*bb;
        x=double(c1+sqrt(aa)*cos(theta_rad));
        y=double(c2+sqrt(bb)*sin(theta_rad));
        [x_cw, y_cw] = poly2cw(x,y);
        [stabset_x, stabset_y] = polybool('subtraction', stabset_x, stabset_y, x_cw, y_cw);
    end
end

sols_output_20=[stabset_x', stabset_y'];
clearvars stabset_x stabset_y;

rerury_matrix = NaN(3,size(sols_output_20,1));
for idx=1:size(sols_output_20,1)
    Ki = sols_output_20(idx,1);
    Kp = sols_output_20(idx,2);
    C = pid(Kp,Ki);
    sys_error = minreal(1/(1+P_tf*C));
    sys_control = minreal(C/(1+P_tf*C));
    sys_output = minreal(P_tf*C/(1 + P_tf*C));
    [re, f_re] = getPeakGain(sys_error);
    [ru, f_ru] = getPeakGain(sys_control);
    [ry, f_ry] = getPeakGain(sys_output);
    rerury_matrix(:,idx) = [re; ru; ry];
end
plot(rerury_matrix(1,1:end-1),rerury_matrix(2,1:end-1));
%plot3(x,y,reru_matrix(3,:),'k-x');
toc