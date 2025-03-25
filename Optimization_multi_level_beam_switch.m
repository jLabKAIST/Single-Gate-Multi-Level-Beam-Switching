clear;close all;clc;
warning('off','all');
addpath("RETICOLO V9/V9/reticolo_allege_v9");
%% Optimization parameters
folderindex = 8000; % save folder index
OptParm.wavelength = 8.0e-6; % wavelength used in geometrical parameter calculation
OptParm.sim_wavelength = 8.0e-6; % wavelength of input wave (useful for broadband spectrum analysis)
OptParm.N = 20; % # of design variables
% For two-level beam switching, set OptParm.EFs = [0.2, 0.6]
% For three-level beam switching, set OptParm.EFs = [0.2, 0.6, 1.0]
% For four-level beam switching, set OptParm.EFs = [0.05, 0.35, 0.65, 0.95]
OptParm.EFs = [0.2, 0.6, 1.0]; % Target Fermi levels
% For two- and three-level beam switching, 
% set OptParm.diffraction_channels = [1 0 -1];
% For four-level beam switching, 
% set OptParm.diffraction_channels = [2 1 0 -1 -2];
OptParm.diffraction_channels = [1 0 -1]; % Propagating diffraction channels
OptParm.gradient_type = 'shape-derivative'; % Either 'grayscale' or 'shape-derivative'
OptParm.angle = 80; % Target diffraction angle
% heights = linspace(1.5*pi,3.5*pi,9); % height sweep
heights = 1.5*pi; 
% thicknesses = linspace(1.5*pi,3.5*pi,9); % thickness sweep
thicknesses = 1.5*pi;
repeats = 1:30;
OptParm.aspect_ratio = 10; % maximum aspect ratio
OptParm.spacer_thickness = 30e-9; % HfO2 spacer thickness
OptParm.t_HfO2 = 50e-9; % HfO2 etch-stop layer thickness
OptParm.b_coefficient = 0.3; % trade-off coefficient
%%

for repeat = repeats
for height_L=heights
   for thickness_L=thicknesses
       clc;
       
       OptParm.height_L = height_L;
       OptParm.thickness_L = thickness_L;
       
       foldername = sprintf('%d_level_gradient_type_%s_angle_%d_height_thickness_sweep_%d_b_coefficient_%.1f', length(OptParm.EFs), OptParm.gradient_type, OptParm.angle, folderindex, OptParm.b_coefficient);
       filename = sprintf('%s/height_%.2f_thickness_%.2f_%d.mat',foldername, height_L,thickness_L, repeat);
       
       if isfolder(foldername) ~= 1
          
           mkdir(foldername);
       end
       
       if isfile (filename) ~= 1
       folderindex
    height_L
    thickness_L
    % clear;


% clear;
close all;

rng('default');
rng('shuffle');
img = rand(1,OptParm.N);
% load('OptimizedStructure\two_level_HfO2_50nm_current.mat','x');
% img = x;
% img=linspace(0,1,N);
lb = zeros(1,OptParm.N);
ub = ones(1,OptParm.N);




objf = @(x) getobjf(x, OptParm);
% global xhistory;
% global fhistory;
% xhistory = [];
% fhistory = [];
% options = optimoptions('fminimax', 'SpecifyObjectiveGradient', true, 'Diagnostics', 'on', 'Display', 'iter', 'OutputFcn', @myoutput);
options = optimoptions('fminimax', 'SpecifyObjectiveGradient', true, 'Diagnostics', 'on', 'Display', 'iter','MaxIterations', 1000);
if strcmp(OptParm.gradient_type,'grayscale')
[x, fval, maxfval, exitflag, output, lambda] = fminimax(objf, img, [], [], [], [], lb, ub, [], options);
else
    %% get feasible starting point
[A, b] = min_feature_size_constraint(OptParm);
proceed_optimization = true;
    for kkk = 1:10000
        if ~all(A*img'-b <=0)
            % disp(img);
            img = rand(1,OptParm.N);
        else
            % disp(kkk);
            proceed_optimization = true;
            break
        end
    end
    %%
    % disp(proceed_optimization);
if ~proceed_optimization
    options_NO = optimoptions('fminimax', 'SpecifyObjectiveGradient', true, 'Diagnostics', 'on', 'Display', 'iter', ...
       'MaxIterations', 1);
[x, fval, maxfval, exitflag, output, lambda] = fminimax(objf, img, A, b, [], [], lb, ub, [], options_NO);
else
[x, fval, maxfval, exitflag, output, lambda] = fminimax(objf, img, A, b, [], [], lb, ub, [], options);
% [x, fval, maxfval, exitflag, output, lambda] = fmincon(objf, img, A, b, [], [], lb, ub, [], options);

end

end

% xsave = x;
% regexprep(num2str(xsave),'\s+',',')


 [T_final, ~] = multi_level_beam_switch_objf(x, OptParm);
if size(OptParm.EFs,2) == 2 %% two-level
 D_final = [T_final(1,1)/sum(T_final(1,:)) T_final(2,3)/sum(T_final(2,:))]; 
elseif size(OptParm.EFs,2) == 3 %% three-level
 D_final = [T_final(1,1)/sum(T_final(1,:)) T_final(2,2)/sum(T_final(2,:)) T_final(3,3)/sum(T_final(3,:))]; 
elseif size(OptParm.EFs,2) == 4 %% four-level
 D_final = [T_final(1,1)/sum(T_final(1,:)) T_final(2,2)/sum(T_final(2,:)) T_final(3,4)/sum(T_final(3,:)) T_final(4,5)/sum(T_final(4,:))]; 
end

 save(filename);
   end
   end
end
end



function [F, G] = getobjf(x, OptParm)

[T, grad_adjs] = multi_level_beam_switch_objf(x, OptParm);
[F, G] = cal_FoM_multi_level(T, grad_adjs, OptParm);

end

function [A, b] = min_feature_size_constraint(OptParm)
if size(OptParm.EFs,2) == 2 %% two-level
period = abs(OptParm.wavelength/sind(OptParm.angle));
elseif size(OptParm.EFs,2) == 3 %% three-level
period = abs(OptParm.wavelength/sind(OptParm.angle));
elseif size(OptParm.EFs,2) == 4 %% four-level
period = 2*abs(OptParm.wavelength/sind(OptParm.angle));
end

OptParm.n_Si = 3.42;
N = OptParm.N;

k = 2*pi/OptParm.wavelength * OptParm.n_Si;
real_height = OptParm.height_L/k;
min_feature_size = real_height/OptParm.aspect_ratio;

dx = period/N*2;
const=min_feature_size*2/dx;

A1 = zeros([N/2 N]);

for i=1:N/2
A1(i,2*i-1) = -1;
A1(i,2*i) = -1;
end

A2 = zeros([N/2 N]);
for i=1:N/2-1
A2(i,2*i)=1;
A2(i,2*i+1)=1;
end
A2(N/2,1) = 1;
A2(N/2,N) = 1;

b1 = ones([N/2 1])*(-const);
b2 = 2 - ones([N/2 1])*const;

A = [A1 ; A2];
b = [b1; b2];




end
% 
% function stop = myoutput(x,optimvalues,state)
% global xhistory;
% global fhistory;
%         stop = false;
%         if isequal(state,'iter')
%           xhistory = [xhistory; x];
%           fhistory = [fhistory; optimvalues.fval'];
%         end
%     end
%     