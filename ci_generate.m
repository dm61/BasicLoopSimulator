function ci = ci_generate(ci_meal) % [(mg/dL)/(5 min)]
% generate ci (carb impact, i.e. insulin counteraction) structure from ci_meal entry

% ci_meal example
% ci_meal.carbs = 100; % total meal carbs [g]
% ci_meal.time = [0 45 90 200 300]'; % meal ci time points [minutes]
% ci_meal.value = [0 1 0 0 0]'; % meal ci relative values, scaled in ci_generate func 
% ci_meal.start = meal_start/5+1; % start time index

% global variables defined in BasicLoopSimulator.m
global ISF; % insulin sensitivity factor [(mg/dL)/U]
global CIR; % carb to insulin ratio [g/U]
global n_sim; % number of 5min steps in the simulation

% setup ci (carb impact) structure
ci.time = (0:5:(n_sim-1)*5)';
ci.value = zeros(1,n_sim)';
    
% aux variables
n_meal_end = max(ci_meal.time)/5+1;
ci_temp.time = (0:5:(n_meal_end-1)*5)';
ci_temp.value = interp1q(ci_meal.time,ci_meal.value,ci_temp.time);
scale = ci_meal.carbs*(ISF/CIR)/sum(ci_temp.value);
ci_temp.value = ci_temp.value.*scale;
       
% calculate carb impact (counteraction) values
if(ci_meal.start+n_meal_end-1 < n_sim)
    ci.value(ci_meal.start:ci_meal.start+n_meal_end-1) = ci_temp.value;
else
    ci.value(ci_meal.start:n_sim) = ci_temp.value(1:n_sim-ci_meal.start+1);
end

end