%% Basic Loop Simulator
% dm61 8/5/2017
% Assumptions: 
%   IOB(0)=0, COB(0)=0
%   Single meal with known carb absorption  
%   Ideal system model: deltaBG = (carb impact)-(insulin impact)
%   Does not incude Loop RC or momentum effects
%   Includes exponential insulin absorption curves with td, tp parameters
%   Includes current Loop and dynamic dosing algorithms, see: algorithm below 
global DIA; % duration of insulin action [h], td = DIA
global ISF; % insulin sensitivity factor [(mg/dL)/U]
global CIR; % carb to insulin ratio [g/U]
global n_sim; % number of 5-min steps
% setup meal, system and simulation parameters
meal_carbs = 50; % grams of carbs in the meal
meal_absorption_time = 6*60; % carbs absorption time in minutes
% default constant absorption rate, arbitrary curve can be setup below
meal_start_time = 60; % meal time in minutes after start of simulation
pre_bolus_time = 20; % pre-bolus time (i.e. time when meal entered), in minutes ahead of the meal
bg_initial = 100; % initial bg value [mg/dL]
% algorithm options
algorithm.bolus = true; % true = use dynamic dosing for bolus, false = current Loop
algorithm.temp = true; % true = use dynamic dosing for temps, false = current Loop
algorithm.accept = true; % true = accept and deliver bolus, false = let Loop handle all
algorithm.alpha = 1; % agressiveness factor, 0 = no dynamic super bolusing
sim_time = 10; % simulation time [h]
pause_times = []; % array of time indeces (1 to n_sim), e.g. [1 10 20] to pause and show current display prediction curve
% system parameters
DIA = 6;
ISF = 50;
CIR = 10;
bg_target_max = 100; % uppper limit of the bg target
bg_target_min = 100; % lower limit of the bg target
bg_guard = 70; % bg guard, no bolus suggested below this level, zero temp below this level
bg_target = (bg_target_max+bg_target_min)/2; % target bg value
basal_rate = 0.5; % nominal basal rate [U/h]
max_temp_rate = 5.0; % maximum temp [U/h]
min_temp_rate = -basal_rate; % minimum temp (negative) [U/h]
% simulation setup
n_sim = round(sim_time*60/5)+1; % total number of simulation points
sim_time_minutes = n_sim*5;
times = 0:5:(n_sim-1)*5;
n_bolus = round((meal_start_time - pre_bolus_time)/5)+1; % time index when meal entered and bolus suggested
nDIA = round(DIA*60/5)+1; % number of time slots in DIA
%% setup insulin absorption curves, bg array, meal 
% normalized scalable exponential insulin activity and IOB curves
td = DIA*60; % insulin duration in minutes
tp = 55; % insulin peak time, nominally Novolog = 75, fiasp = 55
% insulin curves
tau = tp*(td-tp)/(td-2*tp); % time constant of exp decay
S = ((td/tau)^2)*exp(td/tau)/((td-2*tau)*exp(td/tau)+td+2*tau); % aux scale factor
Ia = @(t) S.*(t./td).*(1-t./td).*exp(-t/tau); % insulin activity (AUC=1)
IOB = @(t) 1-(S/td).*(tau/td).*(exp(-t/tau).*(t.^2 - td*(t+tau) +2*tau^2 + 2*tau*t)+tau*(td-2*tau));
%% dynamic dosing algorithm parameters and functions
low_temp_scale = -(min_temp_rate/60)*ISF; % asymptote of the bg rate of change due to zero temping
% bg impact of zero temping as a function of time
BGTempImpact = @(t) low_temp_scale*(S/td)*(tau/td)*tau*...
    (exp(-t/tau).*(-6*tau^2+2*tau*(td-2*t)+t.*(td-t))+...
    6*tau^2-2*tau*td-2*tau*t+td*t);
% dynamic target
BGtarget = @(t) bg_guard + (bg_target-bg_guard)*t/td; % dynamic target
%% initiate bg array
bg = zeros(n_sim,1); % initiate bg array
bg(1) = bg_initial; % initial bg value
bg_predicted = bg; % predicted bg array
%% meal setup
ci_meal.carbs = meal_carbs; % total meal carbs [g]
ci_meal.time = [0 1 meal_absorption_time-1 meal_absorption_time]'; % meal ci time points [minutes]
ci_meal.value = [0 1 1 0]'; % meal ci relative values, scaled in ci_generate func 
meal_start = meal_start_time; % meal start time [min], default to 60min
ci_meal.start = meal_start/5+1; % start time index
ci = ci_generate(ci_meal); % generate actual carb impact array [(mg/dL)/5min]
% what Loop thinks the meal is (assumed ideally the same as of now)
ci_meal_loop = ci_meal; % ideally the same as the actual carbs
meal_start_loop = meal_start; % meal start time (Loop) [min]
ci_meal_loop.start = meal_start_loop/5+1; % (Loop) start time index
ci_loop = ci_generate(ci_meal_loop); % generate Loop model carb impact array
% plot(ci.time,ci.value); check carb impact vector
% grid on;    
%% setup insulin impact array and temp basal array
bg_impact = zeros(n_sim,1);
temp_basal = zeros(n_sim,1);
DIA_times = 5:5:(nDIA-1)*5; % array of DIA times, excluding zero
I_activity_array = 5*Ia(0:5:(nDIA-1)*5)'; % array of insulin activity over DIA
%
%% run simulation 
for i=2:n_sim
    if(i == n_bolus) % enter meal and deliver bolus
        for j=2:n_sim % enter carb impact to bg_predicted array
            bg_predicted(j) = bg_predicted(j-1) + ci.value(j); 
        end
        if(~algorithm.bolus)
            % conventional dosing
            bolus = (bg_predicted(i+nDIA-1)-bg_target)/ISF;
            bolus = max(0,bolus);
        else
            % dynamic dosing
            suggested_bolus_array = ...
                (bg_predicted(i+1:i+nDIA-1)+algorithm.alpha*BGTempImpact(DIA_times)'-BGtarget(DIA_times)')./ ...
                (1-IOB(DIA_times)')/ISF;
            bolus = min(suggested_bolus_array);
            bolus = max(0,bolus);
        end
        if(~algorithm.accept) bolus = 0; end % option to not deliver suggested bolus
        bg_impact(i:i+nDIA-1) = bg_impact(i:i+nDIA-1)+bolus*ISF*I_activity_array;
    end
    bg(i) = bg(i-1) + ci.value(i)- bg_impact(i); % actual bg update
    if( i>= n_bolus)
        bg_predicted = bg; % reset predicted bg values to current bg array
    else
        bg_predicted = bg_initial*ones(n_sim,1);
    end
    if(i+nDIA <= n_sim)
        if(i >= n_bolus) % update prediction only after the meal is entered
            for j=i+1:n_sim % update prediction based on insulin delivered so far
                bg_predicted(j) = bg_predicted(j-1) + ci.value(j)-bg_impact(j); % model based simulation
            end
        end
        if any(pause_times == i)
            clf;
            plot(ci.time,bg_predicted,'k','LineWidth',2);
            disp(bg_predicted(i+nDIA));
            grid on;
            pause;
        end
        if(i >= n_bolus) % temps only if meal has been entered since we assume IOB=0,COB=0 initially 
            if(~algorithm.temp)
                % current Loop temp dosing
                five_minute_dose = ((bg_predicted(i+nDIA)-bg_target)/ISF)*5/30;
            else
                % dynamic temp dosing
                suggested_dose_array = ...
                    (bg_predicted(i+1:i+nDIA-1)-BGtarget(DIA_times)')./ ...
                    (1-IOB(DIA_times)')/ISF;
                suggested_dose = min(suggested_dose_array);
                five_minute_dose = suggested_dose*5/30; % spread over 30 minutes as Loop currently does
                %five_minute_dose = suggested dose; % all in 5 min
            end
        else
            five_minute_dose = 0;
        end
        if(five_minute_dose*60/5 > max_temp_rate)
            five_minute_dose = max_temp_rate*5/60;
        end
        if(five_minute_dose*60/5 < min_temp_rate)
            five_minute_dose = min_temp_rate*5/60;
        end
        if any(bg_predicted(i+1:n_sim) < bg_guard)
            five_minute_dose = min_temp_rate*5/60;
        end
        % add temp effect to insulin impact array
        bg_impact(i+1:i+nDIA) = bg_impact(i+1:i+nDIA)+ ...
            five_minute_dose*ISF*I_activity_array; 
        % update temp basal array (for plotting)
        temp_basal(i+1) = five_minute_dose*60/5; % temp basal [U/h]
    end
end
[xt,yt] = stairs(ci.time,temp_basal); % nicer plot of temps
%
% outputs
fprintf('\ndosing algorithm options\n');
disp(algorithm);
fprintf('bolus = %4.2f\n',bolus); % maximum bg over simulation time
fprintf('maximum BG = %3.0f\n',max(bg)); % maximum bg over simulation time
fprintf('minimum BG = %3.0f\n',min(bg)); % minimum bg over simulation time 
% plot results
figure(1);
clf;
%
subplot(3,1,1) % insulin activity and carb counteraction  
    hold on;
    plot(ci.time-meal_start_time,bg_impact,'g','LineWidth',2);
    plot(ci.time-meal_start_time,ci.value,'r','LineWidth',2);
     axis([-meal_start_time ...
     sim_time_minutes-meal_start_time ...
     -1 inf]);
    grid on;
    hold off;
%
subplot(3,1,2) % bg
    hold on;
    plot(ci.time-meal_start_time,bg,'b','LineWidth',2);
    plot(ci.time-meal_start_time,bg_guard*ones(n_sim,1),'r','LineWidth',1);
    plot(ci.time-meal_start_time,bg_target*ones(n_sim,1),'--k','LineWidth',1);
    grid on;
     axis([-meal_start_time ...
     sim_time_minutes-meal_start_time ...
     50 inf]);
    hold off;
%
subplot(3,1,3) % temps
    hold on;
    plot(xt-meal_start_time,yt,'k','LineWidth',2);
    grid on;
     axis([-meal_start_time ...
     sim_time_minutes-meal_start_time ...
     min_temp_rate max_temp_rate]);
    hold off;