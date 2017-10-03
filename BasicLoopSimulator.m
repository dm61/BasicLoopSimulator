% A Very Basic Loop Simulator
% @dm61 8/6/2017
% 9/3/2017 make code compatible with Octave 

% Assumptions and limitations 
%   Initially IOB(0)=0, COB(0)=0
%   Simulation of a single meal with known carb absorption  
%   Ideal system model: deltaBG = (carb impact)-(insulin impact)
%   Does not incude Loop RC or momentum effects or dynamic carb algorithm
%   Assumes a single fixed nominal basal rate during entire simulation
%   Includes exponential insulin absorption curves with td, tp parameters
%   Includes Loop v1.4 and Loop v1.5 dosing algorithms (see: algorithm section below) 

% global variables accessed from ci_generate
global DIA; % duration of insulin action [h], td = DIA
global ISF; % insulin sensitivity factor [(mg/dL)/U]
global CIR; % carb to insulin ratio [g/U]
global n_sim; % number of 5-min simulation steps

% meal parameters
meal_carbs = 50; % grams of carbs in the meal
CA = 3; % carb absorption time in hours
% default constant absorption rate, arbitrary curve can be setup below
meal_start_time = 60; % meal time in minutes after start of simulation
pre_bolus_time = 20; % pre-bolus time (i.e. time when meal entered), in minutes ahead of the meal
bg_initial = 100; % initial bg value [mg/dL]

% algorithm options
algorithm.bolus = true; % true = use Loop v1.5 dosing for bolus, false = Loop v1.4 dosing
algorithm.temp = true; % true = use Loop v1.5 dosing for temps, false = Loop v1.4 dosing
algorithm.accept = true; % true = accept and deliver bolus, false = let Loop handle all
algorithm.alpha = 0.0; % dynamic super bolusing agressiveness (alpha = 0 to 1), 0 = no dynamic super bolusing (as in Loop v1.5)
algorithm.target = 'dynamic'; % 'dynamic' (blend as in Loop v1.5) or 'target' (least agressive) or 'suspend' (most agressive)
sim_time = 12; % total simulation time [h]
pause_times = []; % array of time indeces (1 to n_sim), e.g. [1 10 20] to pause and display current BG prediction curve

% user selectable model parameters
ISF = 50; % insulin sensitivity [(mg/dL)/U]
CIR = 10; % carb ratio [g/U]
bg_target_max = 100; % uppper limit of the bg correction range
bg_target_min = 100; % lower limit of the bg correction range
bg_suspend = 70; % bg suspend, no bolus suggested below this level, zero temp below this level
bg_target = (bg_target_max+bg_target_min)/2; % target bg value
basal_rate = 0.5; % nominal basal rate [U/h]
max_temp_rate = 10.0; % maximum temp [U/h]

% insulin absorption model
DIA = 6; % DIA = td (it is recommended to keep this as is)
td = DIA*60; % insulin duration in minutes
tp = 75; % insulin peak time, nominally Novolog = 75, FIASP = 55

% simulation setup
meal_absorption_time = CA*60; % carbs absorption time in minutes
min_temp_rate = -basal_rate; % minimum temp (negative) [U/h]
n_sim = round(sim_time*60/5)+1; % total number of simulation points
sim_time_minutes = n_sim*5;
n_bolus = round((meal_start_time - pre_bolus_time)/5)+1; % time index when meal entered and bolus suggested
nDIA = round(DIA*60/5)+1; % number of time slots in DIA

% normalized scalable exponential insulin activity and IOB curves
tau = tp*(td-tp)/(td-2*tp); % time constant of exp decay
S = ((td/tau)^2)*exp(td/tau)/((td-2*tau)*exp(td/tau)+td+2*tau); % aux scale factor
Ia = @(t) S.*(t./td).*(1-t./td).*exp(-t/tau); % insulin activity (AUC=1)
IOB = @(t) 1-(S/td).*(tau/td).*(exp(-t/tau).*(t.^2 - td*(t+tau) +2*tau^2 + 2*tau*t)+tau*(td-2*tau));

% dynamic super bolus (not yet in Loop): zero-temping BG impact functions
low_temp_scale = -(min_temp_rate/60)*ISF; % asymptote of the bg rate of change due to zero temping
% bg impact of zero temping as a function of time
BGTempImpact = @(t) low_temp_scale*(S/td)*(tau/td)*tau*...
    (exp(-t/tau).*(-6*tau^2+2*tau*(td-2*t)+t.*(td-t))+...
    6*tau^2-2*tau*td-2*tau*t+td*t);

% BG target function
% (note: blend is slightly different, needs to be updated to the version in Loop 1.5)
switch algorithm.target
    case 'target' % least agressive
        BGtarget_bolus = @(t) bg_target; % fixed target at bg_target for bolus
        BGtarget_temp = @(t) bg_target; % fixed target at bg_target for temps
    case 'suspend' % most agressive
        BGtarget_bolus = @(t) bg_suspend; % fixed target at bg_suspend for bolus
        BGtarget_temp = @(t) bg_suspend + (bg_target-bg_suspend)*t/(DIA*60); % blend target for temps
    otherwise % dynamic target blend from bg_suspend at t=0 to bg_target at td
        BGtarget_bolus = @(t) bg_suspend + (bg_target-bg_suspend)*t/(DIA*60); % blend target for bolus
        BGtarget_temp = @(t) bg_suspend + (bg_target-bg_suspend)*t/(DIA*60); % blend target for temps
end
        
% meal setup, flat rate assumed below, can be modified 
ci_meal.carbs = meal_carbs; % total meal carbs [g]
ci_meal.time = [0 1 meal_absorption_time-1 meal_absorption_time]'; % meal ci time points [minutes]
ci_meal.value = [0 1 1 0]'; % meal ci relative values, scaled in ci_generate func 
meal_start = meal_start_time; % meal start time [min], default to 60min
ci_meal.start = meal_start/5+1; % start time index
ci = ci_generate(ci_meal); % generate actual carb impact array [(mg/dL)/5min]

% initiate bg array
bg = zeros(n_sim,1); % initiate bg array
bg(1) = bg_initial; % initial bg value
bg_predicted = bg; % predicted bg array

% setup insulin impact array and temp basal array
bg_impact = zeros(n_sim,1);
temp_basal = zeros(n_sim,1);
DIA_times = 5:5:(nDIA-1)*5; % array of DIA times, excluding zero
I_activity_array = 5*Ia(0:5:(nDIA-1)*5)'; % array of insulin activity over DIA

% run simulation 
for i=2:n_sim
    
    % bolus calculation at time slot n_bolus
    if(i == n_bolus) % enter meal and deliver bolus
        for j=2:n_sim % enter carb impact to bg_predicted array
            bg_predicted(j) = bg_predicted(j-1) + ci.value(j); 
        end
        if(~algorithm.bolus)
            % Loop v1.4 dosing
            bolus = (bg_predicted(i+nDIA-1)-bg_target)/ISF;
            bolus = max(0,bolus);
        else
            % Loop v1.5 dosing
            suggested_bolus_array = ...
                (bg_predicted(i+1:i+nDIA-1)+algorithm.alpha*BGTempImpact(DIA_times)'-BGtarget_bolus(DIA_times)')./ ...
                (1-IOB(DIA_times)')/ISF;
            bolus = min(suggested_bolus_array);
            bolus = max(0,bolus);
        end
        if(~algorithm.accept) % algorithm option to not deliver suggested bolus
            bolus = 0; 
        end 
        bg_impact(i:i+nDIA-1) = bg_impact(i:i+nDIA-1)+bolus*ISF*I_activity_array;
    end
    
    % model: update bg 
    bg(i) = bg(i-1) + ci.value(i)- bg_impact(i); % actual bg update

    % initialize predicted bg array
    if( i>= n_bolus)
        bg_predicted = bg; % reset predicted bg values to current bg array
    else
        bg_predicted = bg_initial*ones(n_sim,1);
    end
    
    % temp dosing in time slot i 
    if(i+nDIA <= n_sim) % do only until nDIA before end of simulation

        if(i >= n_bolus) % update prediction only after the meal is entered
            for j=i+1:n_sim % update prediction based on insulin delivered so far
                bg_predicted(j) = bg_predicted(j-1) + ci.value(j)-bg_impact(j); % model based simulation
            end
        end
        
        if any(pause_times == i) % pause to display current prediction, enter to continue
            clf;
            plot(ci.time,bg_predicted,'k','LineWidth',2);
            disp(bg_predicted(i+nDIA));
            grid on;
            pause;
        end
        
        if(i >= n_bolus) % temps only if meal has been entered since we assume IOB=0,COB=0 initially
            if(~algorithm.temp)
                % Loop v1.4 temp dosing
                five_minute_dose = ((bg_predicted(i+nDIA)-bg_target)/ISF)*5/30;
            else
                % Loop v1.5 temp dosing
                suggested_dose_array = ...
                    (bg_predicted(i+1:i+nDIA-1)-BGtarget_temp(DIA_times)')./ ...
                    (1-IOB(DIA_times)')/ISF;
                suggested_dose = min(suggested_dose_array);
                five_minute_dose = suggested_dose*5/30; % spread over 30 minutes as Loop currently does
                %five_minute_dose = suggested dose; % all in 5 min
            end
        else
            five_minute_dose = 0;
        end
        
        % various Loop conditions
        if(five_minute_dose*60/5 > max_temp_rate) % max temp rate
            five_minute_dose = max_temp_rate*5/60;
        end
        
        if(five_minute_dose*60/5 < min_temp_rate) % min temp rate
            five_minute_dose = min_temp_rate*5/60;
        end
        
        if any(bg_predicted(i+1:n_sim) < bg_suspend) % bg suspend
            five_minute_dose = min_temp_rate*5/60;
        end
        
        % add temp effect to insulin impact array
        bg_impact(i+1:i+nDIA) = bg_impact(i+1:i+nDIA)+ ...
            five_minute_dose*ISF*I_activity_array; 
        
        % update temp basal array (for plotting)
        temp_basal(i+1) = five_minute_dose*60/5; % temp basal [U/h]
    
    end
end

% outputs
[BGmin, indexBGmin] = min(bg);
[BGmax, indexBGmax] = max(bg);
fprintf('\ndosing algorithm options\n');
disp(algorithm);
fprintf('bolus = %4.2f\n',bolus); % maximum bg over simulation time
fprintf('maximum BG = %3.0f\n',BGmax); % maximum bg over simulation time
fprintf('minimum BG = %3.0f\n',BGmin); % minimum bg over simulation time 
plot_times = ci.time-meal_start_time;
[xt,yt] = stairs(ci.time,temp_basal); % nicer plot of temps

% various plot labels
strmin = num2str(round(BGmin));
strmax = num2str(round(BGmax));
if(algorithm.bolus)
    titlestr = ['Loop v1.5 dosing, \alpha = ' num2str(algorithm.alpha)];
else
    titlestr = 'Loop v1.4 dosing';
end
if(algorithm.accept)
    bolusstr = [', bolus = ' num2str(round(bolus*100)/100)];
else
    bolusstr = ', skip bolus';
end
if(algorithm.temp)
    tempstr = ['Loop v1.5 dosing, min = ' num2str(min_temp_rate) ', max = ' num2str(max_temp_rate) ];
else
    tempstr = ['Loop v1.4 dosing, min = ' num2str(min_temp_rate) ', max = ' num2str(max_temp_rate) ];
end

% plot results
% choose where to plot depending on algorithm choice
if(algorithm.bolus)
    if(algorithm.temp)
        figure(1);
    else
        figure(2);
    end
else
    figure(3);
end
clf;

subplot(3,1,1) % plot insulin activity and carb counteraction  
    hold on;
    plot(plot_times,bg_impact,'g','LineWidth',2);
    plot(plot_times,ci.value,'r','LineWidth',2);
     axis([-meal_start_time ...
     sim_time_minutes-meal_start_time ...
     -1 inf]);
    title(['Insulin action and counteraction, ' titlestr bolusstr ...
        ', DIA = ' num2str(DIA) ', tp = ' num2str(tp)],...
        'FontSize',10,'FontWeight','normal');
    grid on;
    hold off;

subplot(3,1,2) % bg
    hold on;
    plot(plot_times,bg,'b','LineWidth',2);
    plot(plot_times,bg_suspend*ones(n_sim,1),'r','LineWidth',1);
    plot(plot_times,bg_target*ones(n_sim,1),'--k','LineWidth',1);
    grid on;
     axis([-meal_start_time ...
     sim_time_minutes-meal_start_time ...
     50 200]);
    title('BG [mg/dL]','FontSize',10,'FontWeight','normal');
    text(plot_times(indexBGmin),BGmin,strmin);
    text(plot_times(indexBGmax),BGmax,strmax);
    hold off;

subplot(3,1,3) % temps
    hold on;
    plot(xt-meal_start_time,yt,'k','LineWidth',2);
    grid on;
     axis([-meal_start_time ...
     sim_time_minutes-meal_start_time ...
     min_temp_rate max_temp_rate]);
    title(['Basal [U/h], ' tempstr], ...
        'FontSize',10,'FontWeight','normal');
    hold off;
    