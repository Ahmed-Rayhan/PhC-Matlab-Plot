clc;
clear all;
close all;

n_Si = 3.465;    
n_SiO2 = 1.442;  
n_water = 1.333;          
n_cholera = 1.365;        
n_ecoli = 1.388;         
n_shigella = 1.422;       

lambda_design = 1700; 
d_Si = lambda_design / (4 * n_Si);    
d_SiO2 = lambda_design / (4 * n_SiO2);
d_defect = d_Si + d_SiO2; 

lambda_range = linspace(1100, 2750, 6603); 

N = 10; 

transmittance_water = arrayfun(@(lambda) transfer_matrix_method(lambda, N, n_Si, n_SiO2, n_water, d_Si, d_SiO2, d_defect), lambda_range);
transmittance_cholera = arrayfun(@(lambda) transfer_matrix_method(lambda, N, n_Si, n_SiO2, n_cholera, d_Si, d_SiO2, d_defect), lambda_range);
transmittance_ecoli = arrayfun(@(lambda) transfer_matrix_method(lambda, N, n_Si, n_SiO2, n_ecoli, d_Si, d_SiO2, d_defect), lambda_range);
transmittance_shigella = arrayfun(@(lambda) transfer_matrix_method(lambda, N, n_Si, n_SiO2, n_shigella, d_Si, d_SiO2, d_defect), lambda_range);

figure;
plot(lambda_range, transmittance_water, 'b'); hold on;
plot(lambda_range, transmittance_cholera, 'r');
plot(lambda_range, transmittance_ecoli, 'g');
plot(lambda_range, transmittance_shigella, 'm');
xlabel('Wavelength (nm)');
ylabel('Transmittance');
legend('Pure Water', 'Vibrio cholera', 'Escherichia coli', 'Shigella flexneri');
title('Transmission Spectra for Different Water Samples');
grid on;


figure;
plot(lambda_range, transmittance_water, 'b'); hold on;
plot(lambda_range, transmittance_cholera, 'r');
plot(lambda_range, transmittance_ecoli, 'g');
plot(lambda_range, transmittance_shigella, 'm');
xlim([1800 1900]); 
xlabel('Wavelength (nm)');
ylabel('Transmittance');
legend('Pure Water', 'Vibrio cholera', 'Escherichia coli', 'Shigella flexneri');
title('Enlarged View of Defect Modes');
grid on;

function T = transfer_matrix_method(lambda, N, n_Si, n_SiO2, n_defect, d_Si, d_SiO2, d_defect)
    n_air = 1; 
    n_substrate = 1.517; 

    M_total = eye(2);
    
    for i = 1:N/2
        
        delta_Si = 2 * pi * n_Si * d_Si / lambda;
        M_Si = [cos(delta_Si), 1i * sin(delta_Si) / n_Si; 1i * n_Si * sin(delta_Si), cos(delta_Si)];
        
        delta_SiO2 = 2 * pi * n_SiO2 * d_SiO2 / lambda;
        M_SiO2 = [cos(delta_SiO2), 1i * sin(delta_SiO2) / n_SiO2; 1i * n_SiO2 * sin(delta_SiO2), cos(delta_SiO2)];
        
        M_total = M_Si * M_SiO2 * M_total;
    end
    
    delta_defect = 2 * pi * n_defect * d_defect / lambda;
    M_defect = [cos(delta_defect), 1i * sin(delta_defect) / n_defect; 1i * n_defect * sin(delta_defect), cos(delta_defect)];
    
    M_total = M_total * M_defect * M_total;
    
    A = M_total(1,1);
    B = M_total(1,2);
    C = M_total(2,1);
    D = M_total(2,2);
    
    T = (4 * n_air * n_substrate) / abs(A * n_substrate + B + C * n_air * n_substrate + D * n_air)^2;
end
