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

z=1.5; % defocus in um
sxy=.04095; % um per pixel
fpropmat=exp(kz_of_kxky*(1i*klen/sxy*z));
        
%% look at ft, make sure the fresnel plate is sampled with sufficient
%  resolution so that fresnel transform will work
ift(fpropmat)
%% the 1d otf had values up to pixel 41 of 100
solidpupil=tmp2>0;


%% copy averaged psf into center of a big array
zpos=[-25,-20,-15,0,15,20,25]';
pupils=repmat(tmp2>0,[1,1,length(zpos)])*1i;
avgpsf=newim(s(1),s(2),length(zpos))+.0001;
fpropmat=avgpsf.*1i;
for i=1:length(zpos)
    a=readim(sprintf('/home/martin/1007/avgpsf/psf%+03d.ics',zpos(i)));
    [row,col]=size(a);
    % avgpsf(:,:,i-1)=dip_resampling(a,[s(1)/row s(2)/col],[0 0],'3-cubic');
    avgpsf(512-59:512+60,512-59:512+60,i-1)=DampEdge(sqrt(a(:,:)),.7,2,0);
    fpropmat(:,:,i-1)=exp(kz_of_kxky*(1i*klen/sxy*zpos(i)/10));
end

%%
niter=25;
pupillist=newim(410,410,niter)*1j;
for iter=0:niter-1
    sample=dip_fouriertransform(pupils.*fpropmat,'forward',[1 1 0]);
    sample=sample./abs(sample) .* avgpsf;
    ksample=dip_fouriertransform(sample,'inverse',[1 1 0]);
    ksample=ksample ./ fpropmat;
    current_pupil=dip_mean(ksample,repmat(solidpupil,[1 1 length(zpos)]),[0 0 1]);
    pupils=repmat(current_pupil,[1 1 length(zpos)]);
    pupillist(:,:,iter)=extract(current_pupil,[410 410]);
end
%%



%%


niter=100;
pupillist=newim(410,410,niter)*j;
for i=0:niter-1
    b2=ft(fpropmat.*pupil);
    b2=b2./abs(b2).*b;
    pupil=ift(b2)./fpropmat.*solidpupil;
    pupillist(:,:,i)=extract(pupil,[410 410]);
end

%%
ft(fpropmat.*pupil)