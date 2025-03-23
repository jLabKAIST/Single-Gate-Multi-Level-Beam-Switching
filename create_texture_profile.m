function OptParm = create_texture_profile(OptParm)

OptParm.n_sub = OptParm.n_Si;
OptParm.n_grating = OptParm.n_Si;
OptParm.n_spacer = OptParm.n_SiN;
OptParm.n_BR = OptParm.n_Au;




OptParm.k_pillar = 2*pi/(OptParm.wavelength/real(OptParm.n_grating));
OptParm.k_sub = 2*pi/(OptParm.wavelength/real(OptParm.n_sub));
OptParm.height  = OptParm.height_L/OptParm.k_pillar;
OptParm.thickness =  OptParm.thickness_L/OptParm.k_sub;
if size(OptParm.n_grs,2)==2 %% two-level
OptParm.period = abs(OptParm.wavelength/sind(OptParm.angle));
elseif size(OptParm.n_grs,2)==3 %% three-level
OptParm.period = abs(OptParm.wavelength/sind(OptParm.angle));
elseif size(OptParm.n_grs,2)==4 %% four-level
OptParm.period = 2*abs(OptParm.wavelength/sind(OptParm.angle));
end

% textures for all layers including the top and bottom layers
textures =cell(1,5);
textures{1}= OptParm.n_air; % Uniform, top layer

if strcmp(OptParm.gradient_type,'grayscale')
dx = OptParm.period/OptParm.N;
x = [1:OptParm.N]*dx - 0.5*OptParm.period;
x_center = x - OptParm.period/(2*OptParm.N);
nvec = OptParm.img*(OptParm.n_grating - OptParm.n_air) + OptParm.n_air;
nvec = real(nvec);
% Span for mid-layer patterns is [-period/2, period/2]
textures{2}={x, nvec};
else
    % Span for mid-layer patterns is [-period/2, period/2]
  dx = OptParm.period/OptParm.N*2;
x = [1:OptParm.N/2]*dx - 0.5*OptParm.period;
x_center = x - OptParm.period/(OptParm.N);
x_center = [x_center ; x_center];
x_center = reshape(x_center, 1,[]);
    signs = (-1).^(1:OptParm.N);
    x = x_center+signs*(dx/2).*OptParm.img;

    nvec = ((signs+1)/2)*(OptParm.n_grating-OptParm.n_air)+OptParm.n_air;
    textures{2} = {x, nvec};
    

end

%% graphene 위에 Hafnia 를 추가. SiN을 Hafnia 로 변경 


textures{3} = OptParm.n_HfO2;
textures{4} = OptParm.n_gr;
textures{5} = OptParm.n_HfO2;
textures{6}= OptParm.n_sub;
% textures{5}= 1.5;


% n_Au = 3.7923+54.663i;
textures{7}= OptParm.n_BR; % Uniform, bottom layer00000000000000000000000000000000000000000000000000000
% textures{6}= OptParm.n_sub; % Uniform, bottom layer
% textures{6} =  1.5;
 
% profile = {[2*thickness, thickness, 2*thickness], [1, 2, 3]}; %Thicknesses of the layers, and layers, from top to bottom.
OptParm.profile = {[0, OptParm.height,OptParm.t_HfO2, 0.34e-9, OptParm.spacer_thickness, ...
    OptParm.thickness,  0], [1, 2, 3, 4, 5, 6, 7]}; %Thicknesses of 
% OptParm.profile = {[0, OptParm.height, 0.34e-9, OptParm.spacer_thickness, ...
%     0], [1, 2, 3, 4, 5]}; %Thicknesses of 

OptParm.textures = textures;


 OptParm.nvec = nvec;









end