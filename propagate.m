%% 
% all lengths are in micrometer
NA=1.38; % numerical aperture
ri=1.515;% refractive index
lambda=.5; % emission vacuum wavelength
rbfp=NA/lambda; % radius of bfp for oil immersion objective
s=[1024,1024]; % image size, filling the full diameter of the bfp
apertureX=2*rbfp; % size of the bfp image in um
apertureY=2*rbfp;
cosTheta=newim(s)+1;                  
kz_of_kxky=cosTheta;

klen=ri/NA; % radius of ewald sphere in units of the aperture

kx=xx(s,'freq')*(apertureX);                   % normalized to ellipse
ky=yy(s,'freq')*(apertureY);                   % normalized to ellipse               
rho2=kx^2+ky^2;

tmp2=klen^2-rho2;   % k_0 here is the ratio between k_0 in pixels and aperture in pixels
cosTheta(tmp2>0)=sqrt(tmp2(tmp2>0));        
                    % Theta, being the angle to the optic axis
kz_of_kxky=cosTheta;% still in units of kx and ky

z=3.; % defocus in um
sxy=.04095; % um per pixel
fpropmat=exp(kz_of_kxky*(1i*klen/sxy*z));
        % is klen in units of aperture?
        % klen contains the refractive index of the immersion medium
        
%% look at ft, make sure the fresnel plate is sampled with sufficient
%  resolution so that fresnel transform will work
ift(fpropmat)


%% define a ring in the bfp
pupil=rr(s)>20 & rr(s)<30;
ft(fpropmat.*pupil)