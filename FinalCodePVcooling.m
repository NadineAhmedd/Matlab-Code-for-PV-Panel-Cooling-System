%% constants
irradiance=[770,910,1100,1300,990,870,690];       %irradiance along 7 hours of the day
T_before=[55 58 63 59 58 56 55];                  %measured panel temperature along 7 hours before cooling
time=[9:15];                                      %from 9 AM till 3 PM
T_ref=25;                                         %refrence temperature at STP (C)
Area=1485*668*1e-6;                               %PV panel area (m^2)
%for porous material
L=0.05;                                           %panel height(m)
Keq=0.693;                                        %equivalent thermal conductivity(W/(m*K))
q_dot=200;                                        %heat flux density (W/m^2)
h = 20;                                           %heat transfer coefficient(W.m^-2.K)
T_coef=0.0045;                                    %temperature coefficient
n_ref=16.677;                                     %refrence efficiency(%)
T_amb_values=[36 38 39 40 40 41 41];              %ambient temperature along 7 hours of the day
%for nanofluid
Temp_NF = [35.5593 38.1814 40.59095 43.75 46 37.5 35]; %Temperature values after MgO_T
%for heatpipe
T_coeff=-0.45; 
nominal_power=330;
eff_panel=16.677/100;
Temp_HPT=[50   53   57   55   52   51   54  ];

%% Temperature calculations for porous material
 for i = 1:length(T_amb_values)
    T_amb = T_amb_values(i);
    % Calculate T_0 and T_L using the given equations
    T_0=(h^2*L*T_amb+q_dot*h*L+h*Keq*T_amb+q_dot*Keq)/(2*h*Keq+h^2*L);
    T_L=(q_dot*h*L+h^2*L*T_amb+q_dot*h*L+2*h*Keq*T_amb+q_dot*Keq)/(2*h*Keq+h^2*L);
    % Store values in arrays
    T_0_values_Porous(i) = T_0;
    T_L_values_Porous(i) = T_L;
 end
 
%% power calculations 
%for porous material
power=irradiance*n_ref;
powerloss_before=(T_before-T_ref)*T_coef;
power_before=((1-abs(powerloss_before)).*(power))/100;
powerloss_porous_back=(T_0_values_Porous-T_ref)*T_coef;
PowerPorous=((1-abs(powerloss_porous_back)).*power)/100;
%for nanofluid
% Calculate power output for each irradiance value and temperature
PowerNanoFluid = irradiance .* (effic_r * (1 - 0.0045 * (Temp_NF - Tr_nf)));
power =irradiance*eff_panel ;
Per_of_powerloss_HPT=(Temp_HPT-T_ref)*T_coeff; % with HPT
PowerHPT=(1-abs(Per_of_powerloss_HPT/100)) .* (power);
Per_of_powerloss=(T_before-T_ref)*T_coeff; %without coolling system
power_without=(100-abs(Per_of_powerloss)).* (power/100);
diff=PowerHPT-power_without;

%% Efficiency Calculations
efficiency_before=(power_before./(irradiance.*Area))*100;
%for porous material
efficiency_after_Porous=(PowerPorous./(irradiance.*Area))*100;
%for nanofluid
efficiency_after_nf = (PowerNanoFluid./ (irradiance .* Area)) * 100;
%for heatpipe
efficiency_after_HP=(PowerHPT./(irradiance.*Area))*100;

%% Plotting
figure;
plot(time,T_before,'o-','LineWidth' ,2);
hold on;
plot(time,T_0_values_Porous,'o-','LineWidth' ,2);
plot(time,Temp_NF,'o-','LineWidth' ,2);
plot(time,Temp_HPT,'o-','LineWidth' ,2);
legend('Panel Temprature before cooling','Panel Temperature Using Porous Material','Panel Temperature Using NanoFluid MgO','Panel Temperature Using HeatPipe');
xlabel('Time','fontWeight','bold');
ylabel('Temperature (°C)','fontWeight','bold');
title('Temperature Improvement after using Porous Material','fontWeight','bold');
figure;
plot(time,power_before,'o-','LineWidth' ,2);
hold on;
plot(time,PowerPorous,'o-','LineWidth' ,2);
plot(time, PowerNanoFluid, 'o-', 'LineWidth', 2);
plot(time,abs(PowerHPT),'o-','LineWidth', 2);
xlabel('Time','fontWeight','bold');
ylabel('Output Power (Watt)','fontWeight','bold');
title('Output Power Improvement','fontWeight','bold');
legend('Power Before Cooling','Power After Cooling Using Porous Material','Power After Cooling Using NanoFluid MgO','Power After Cooling Using Heat Pipe');
figure;
plot(time,efficiency_before,'o-','LineWidth' ,2);
hold on;
plot(time,efficiency_after_Porous,'o-','LineWidth' ,2);
plot(time,efficiency_after_nf, 'o-', 'LineWidth', 2);
plot(time,efficiency_after_HP,'o-','LineWidth', 2);
legend('Efficiency Before Cooling','Efficiency After Cooling Using Porous Material','Efficiency After Cooling Using NanoFluid MgO','Efficiency After Cooling Using Heat Pipe');
xlabel('Time','fontWeight','bold');
ylabel('Efficiency %','fontWeight','bold');
title('Efficiency Improvement','fontWeight','bold');
%% Enhancement Calculations
%for porous material
disp(['Average Power Before Cooling=',num2str(mean(power_before))]);
disp(['Average Power After Using Porous Material=',num2str(mean(PowerPorous))]);
disp(['Average Efficiency Before Cooling=',num2str(mean(efficiency_before))]);
disp(['Average Efficiency After Using Porous Material=',num2str(mean(efficiency_after_Porous))]);
disp(['Average Temperature After Using Porous Material=',num2str(mean(T_0_values_Porous))]);
disp(['Average Temperature Before Cooling=',num2str(mean(T_before))]);
disp(['Power Improvement Percentage Using Porous Material=',num2str(((mean(PowerPorous)-mean(power_before))/mean(power_before))*100)]);
disp(['Efficiency Improvement Percentage Using Porous Material=',num2str(((mean(efficiency_after_Porous)-mean(efficiency_before))/mean(efficiency_before))*100)]);
disp(['Temperature Enhancement Using Porous Material=',num2str(mean(T_before)-mean(T_0_values_Porous))]);
disp(['Temperature Enhancement Using Porous Material Percentage=',num2str((mean(T_before)-mean(T_0_values_Porous))/mean(T_before)*100)]);
%for nanofluid
disp(['Average Power After Using Nano Fluid=',num2str(mean(PowerNanoFluid))]);
disp(['Average Efficiency After Using NanoFluid=',num2str(mean(efficiency_after_nf))]);
disp((mean(efficiency_after_nf)-mean(efficiency_before))/mean(efficiency_before));

%% for heatpipe
avr_eff_HPT=mean(efficiency_after_HP)
avr_eff_without=mean(efficiency_before)
enhancement_eff=avr_eff_HPT-avr_eff_without
per_eff_enhancement1=(enhancement_eff/avr_eff_without)*100
avr_eff_power=mean(PowerHPT )
avr_eff_without=mean(power_before)
enhancement_power=avr_eff_power-avr_eff_without
per_eff_enhancement2=(enhancement_power/avr_eff_power)*100
