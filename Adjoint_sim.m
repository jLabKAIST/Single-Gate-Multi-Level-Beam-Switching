function [Ex, Ey, Ez] = Adjoint_sim(OptParm, textures, profile)


N = OptParm.N;
wavelength = OptParm.wavelength;
period = OptParm.period;
theta0 = OptParm.theta0;
ordernumber = OptParm.ordernumber;
wavedirection = OptParm.wavedirection;
gratinglayer = OptParm.gratinglayer;
polarization = OptParm.polarization;
gradient_type = OptParm.gradient_type;
textures = OptParm.textures;
profile = OptParm.profile;
%% wavedirection 1: inc_top_transmitted, 2: inc_top_reflected, 3: inc_bottom_transmitted, 4: inc_bottom_reflected




nn = OptParm.nn;
parm = res0(polarization);
if wavedirection==1 || wavedirection == 4
    n_inc = textures{end};
else
    n_inc = textures{1};
end
k_parallel = n_inc*sin(theta0/180*pi);


parm.res1.champ = 1;



% parm.res1.trace = 1;
aa = res1(wavelength, period, textures, nn, k_parallel, parm);
% k_parallel
result = res2(aa, profile);

% orders = result.inc_bottom_transmitted.order
% 
% K = result.inc_bottom_transmitted.K;
% kx = K(:,1);
% angle = asin(kx)/pi*180
    


% parm.res3.trace = 1;
% parm.res3.npts = [101,101,101];



if wavedirection==1


parm.res3.sens = -1;

if polarization == 1
inc = result.inc_bottom.PlaneWave_E(2);
else 
    inc = result.inc_bottom.PlaneWave_H(2);
end
elseif wavedirection==2

parm.res3.sens = 1;
if polarization == 1
inc = result.inc_top.PlaneWave_E(2);
else 
    inc = result.inc_top.PlaneWave_H(2);
end
elseif wavedirection==3

parm.res3.sens = 1;
if polarization == 1
inc = result.inc_top.PlaneWave_E(2);
else 
    inc = result.inc_top.PlaneWave_H(2);
end
else

parm.res3.sens = -1;
if polarization == 1
inc = result.inc_bottom.PlaneWave_E(2);
else 
    inc = result.inc_bottom.PlaneWave_H(2);
end
end


parm.res3.npts = zeros(length(textures),1);
parm.res3.npts(gratinglayer) = 51;
% parm.res3.npts = [51,51,51];

%% gradient_type : 'grayscale' 일 경우
if strcmp(OptParm.gradient_type,'grayscale')
x_dc = cell2mat(textures{gratinglayer}(1));
x = x_dc - period / (2*N);
Nx = 21;
fine_x = [];
temp = linspace(-period/2, x_dc(1),Nx+2);
fine_x = [fine_x temp(2:end-1)];
for i = 1: N-1
temp = linspace(x_dc(i), x_dc(i+1), Nx+2);
fine_x = [fine_x temp(2:end-1)];
end
else
    fine_x = cell2mat(textures{gratinglayer}(1));
end
% x = x_dc;
% inc = result.inc_bottom.PlaneWave_H(2);

% inc = result.inc_bottom.H(2)
% inc = 1;
% x = linspace(-period, period, 101);
% parm.res3.trace=1;
% [e, z, index, wZ,loss_per_layer,loss_of_Z,loss_of_Z_X,X,wX,Flux_Poynting] = res3(fine_x,aa,profile,inc,parm);
[e,z,index] = res3(fine_x,aa,profile,inc,parm);
if polarization == -1
Ex = e(:,:,2);
Ey = [];
Ez = e(:,:,3);
else
    Ex = [];
    Ey = e(:,:,1);
    Ez = [];
    
end
% Ey = mean(e(2:end-1,:,1),1);
% Ey = reshape(Ey, Nx, []);
% Ey = mean(Ey, 1);
% e = e.*sqrt(index);
% Ex = mean(e(2:end-1,:,2),1) ;
% Ez = mean(e(2:end-1,:,3),1);
% Ex = reshape(Ex, Nx,[]);
% Ex = mean(Ex,1);
% Ez = reshape(Ez, Nx, []);
% Ez = mean(Ez,1);

retio;
end


