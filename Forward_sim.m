function [t, Ex, Ey, Ez, angle, Eff, varargout] = Forward_sim(OptParm)

nout = max(nargout,1)-6;

N = OptParm.N;
wavelength = OptParm.wavelength;
sim_wavelength = OptParm.sim_wavelength;
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
%% gradient_type 'shape': shape derivative, 'grayscale': grayscale derivative                                      


nn = OptParm.nn;
parm = res0(polarization);
if wavedirection==1 || wavedirection==2
n_inc = textures{1};
else
    
    n_inc = textures{end};
end
k_parallel = n_inc*sin(theta0/180*pi);

parm.res1.champ = 1;



% parm.res1.trace = 1;
if nout == 0
aa = res1(sim_wavelength, period, textures, nn, k_parallel, parm);
elseif nout ==1
[aa, n_eff] = res1(sim_wavelength, period, textures, nn, k_parallel, parm);
varargout = n_eff{gratinglayer};
end

result = res2(aa, profile);
% disp(sprintf("Number of orders: %d", length(result.inc_top_reflected.efficiency)));
if wavedirection==1
orders = result.inc_top_transmitted.order;
elseif wavedirection==2
    orders = result.inc_top_reflected.order;
elseif wavedirection==3
    orders = result.inc_bottom_transmitted.order;
else
     orders = result.inc_bottom_reflected.order;
end


if isempty(find(orders==ordernumber))
    error('Ordernumber does not exist')
else
    orderindex = find(orders==ordernumber);
    if wavedirection==1
    K = result.inc_top_transmitted.K;
elseif wavedirection==2
    K = result.inc_top_reflected.K;
elseif wavedirection==3
    K = result.inc_bottom_transmitted.K;
else
       K = result.inc_bottom_reflected.K;
end

%     asin(K(:,1))/pi*180;
    kx = K(orderindex,1);
    angle = asin(kx)/pi*180;
    
end
% kx


if wavedirection==1
t = result.inc_top_transmitted.amplitude{ordernumber};
Eff = result.inc_top_transmitted.efficiency{ordernumber};

parm.res3.sens = 1;
if polarization == 1
inc = result.inc_top.PlaneWave_E(2);
else 
    inc = result.inc_top.PlaneWave_H(2);
end

elseif wavedirection==2
t = result.inc_top_reflected.amplitude{ordernumber};
Eff = result.inc_top_reflected.efficiency{ordernumber};
parm.res3.sens = 1;
if polarization == 1
inc = result.inc_top.PlaneWave_E(2);
else 
    inc = result.inc_top.PlaneWave_H(2);
end
elseif wavedirection==3
   t = result.inc_bottom_transmitted.amplitude{ordernumber};
Eff = result.inc_bottom_transmitted.efficiency{ordernumber};
parm.res3.sens = -1;
if polarization == 1
inc = result.inc_bottom.PlaneWave_E(2);
else 
    inc = result.inc_bottom.PlaneWave_H(2);
end
else
t = result.inc_bottom_reflected.amplitude{ordernumber};
Eff = result.inc_bottom_reflected.efficiency{ordernumber};
parm.res3.sens = -1;
if polarization == 1
inc = result.inc_bottom.PlaneWave_E(2);
else 
    inc = result.inc_bottom.PlaneWave_H(2);
end
end





parm.res3.npts = zeros(length(textures),1);
parm.res3.npts(gratinglayer) = 51;



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

% inc = result.inc_bottom.H(2)
% inc = 1;
% x = linspace(-period, period, 101);
% parm.res3.trace=1;
% [e, z, index, wZ,loss_per_layer,loss_of_Z,loss_of_Z_X,X,wX,Flux_Poynting] = res3(fine_x,aa,profile,inc,parm);


% % 시각화용 섹션 - 미사용 시 주석처리
% parm.res3.trace=1;
% parm.res3.npts = [201,51,51,51,51];
% fine_x = linspace(-period/2,period/2,401);
% profile{1}(1) = 10e-6;
% profile{1}(end) = 0e-6;
% %

[e,z,index] = res3(fine_x,aa,profile,inc,parm);
% Forward_sim_draw_Efield(e,z,fine_x,profile,index);
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
% Ex = mean(e(2:end-1,:,2),1);
% Ez = mean(e(2:end-1,:,3),1);
% Ex = reshape(Ex, Nx,[]);
% Ex = mean(Ex,1);
% Ez = reshape(Ez, Nx, []);
% Ez = mean(Ez,1);

% pcolor(fine_x, 

retio;
end


