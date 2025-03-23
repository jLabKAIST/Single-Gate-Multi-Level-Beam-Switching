function OptParm = construct_OptParm(OptParm)
%% Physical and structural parameters
c = 299792458;
OptParm.theta0 = 0;

OptParm.polarization = 1;


OptParm.wavedirection = 2;
 OptParm.gratinglayer = 2;
  OptParm.nn = 40;
%% Material parameters
%%%%%%%%%%%%  Gold -> Babar Gold

Au_nk = table2array(readtable('Material_Index\Au_nk_Babar.csv'));
wvls = Au_nk(:,1)*1e-6;
Au_eps = (Au_nk(:,2)+1j*Au_nk(:,3)).^2;
% Au_eps = real(Au_eps);
Au_nk = sqrt(Au_eps);
OptParm.n_Au = interp1(wvls, Au_nk, OptParm.sim_wavelength);
% disp(OptParm.n_Au^2);
% figure();
% hold on;
% plot(wvls, real(Au_nk));
% plot(wvls, imag(Au_nk));
%%%%%%%%%%%%%
%%%%%%%%%%%% Silver
Ag_nk = table2array(readtable('Material_Index\Ag_nk_Babar.csv'));
wvls = Ag_nk(:,1)*1e-6;
Ag_eps = (Ag_nk(:,2)+1j*Ag_nk(:,3)).^2;
% Au_eps = real(Au_eps);
Ag_nk = sqrt(Ag_eps);
OptParm.n_Ag = interp1(wvls, Ag_nk, OptParm.sim_wavelength);
% figure();
% hold on;
% plot(wvls, real(Ag_nk));
% plot(wvls, imag(Ag_nk));
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% PEC
% OptParm.n_PEC = inf;
%%%%%%%%%%%%%%


OptParm.n_Si = 3.42;
OptParm.n_air = 1;

re_eps = load('Material_Index\[Si3N4_LPCVD2] re_eps.txt');
im_eps = load('Material_Index\[Si3N4_LPCVD2] im_eps.txt');
SiN_nk = sqrt(re_eps(:,2)+1j*im_eps(:,2));
freqs = re_eps(:,1)*1e12;
wvls = c./freqs;
OptParm.n_SiN = interp1(wvls, SiN_nk, OptParm.sim_wavelength);
OptParm.n_SiN_at_target_wvl = interp1(wvls, SiN_nk, OptParm.wavelength);
% disp(OptParm.n_SiN);
OptParm.n_HfO2 = 1.81; % n_HfO2 @ 8 um
% OptParm.n_HfO2 = 1.876; % n_HfO2 @ 7 um
% OptParm.n_HfO2 = 1.8446; % n_HfO2 @ 7.5 um
% OptParm.n_HfO2 = 1.769; % n_HfO2 @ 8.5 um
% OptParm.n_HfO2 = 1.7232; % n_HfO2 @ 9 um


% figure();
% hold on;
% plot(wvls, real(SiN_nk));
% plot(wvls, imag(SiN_nk));

wholeepsilons = load('Material_Index\wholeepsilons.txt');
for i = 1:length(OptParm.EFs)
OptParm.n_grs(i) = sqrt(griddata(wholeepsilons(:,1), wholeepsilons(:,2), wholeepsilons(:,3)+1j*wholeepsilons(:,4), c/OptParm.sim_wavelength,OptParm.EFs(i)));
end

end