% coherent image of mma after fourier filter

%%
s=[64 64]; % calculation size (mma and padding)
zoom=5; % pixels per micro mirror
pad=4; % 4 micromirrors padding
% tm contains the image that is displayed on the mma, tm=0 .. mirror is
% titlted to blaze angle, tm=1 .. mirror is flat

% tm=(rr(s)<4.5); %.*mod(yy(s),2) % hard edged circle
tm=exp(-((xx(s)-14).^2+yy(s).^2)/30); % soft gaussian
% tm=newim(s); % blank
dir=2*(mod(yy(s),2)-.5); % mirror rows tilt in different directions
bigtilt=incimate((1-tm).*dir,zoom);


% for odd zoom, left edge of each mirror is -.25 and right edge +.25
shif=floor(zoom/2); 
if mod(zoom,2)==0
    shif=shif-1;
end
maxtilt=(mod(xx(bigtilt,'corner'),zoom)-shif)/(floor(zoom/2)*4);

mma_phase=newim(maxtilt).*1i;
p=zoom*(pad-1);
q=zoom*(s(2)+1-pad)-1;
% for the simulation the pathdifference should be 2pi for blazing
mma_phase(p:q,p:q)=exp(1i*4*pi*...
    bigtilt(p:q,p:q).*maxtilt(p:q,p:q));

writeim(255*normalize(phase(mma_phase)),'00mma_phase.jpg','JPEG');

%% prepare illumination angles as defined by circular aperture in plane 0
ill_s=[13 13];
ill_kx=xx(ill_s,'freq')/(2*zoom);
ill_ky=yy(ill_s,'freq')/(2*zoom);
% only images inside this circle will contribute to incoherent summation
ill_on=rr(ill_s,'freq')<=.8*.5;
%% for each illumination angle store the intensity in planes 2 and 3
vol_i3=newim([size(maxtilt) ill_s],'double');
vol_i2=newim([size(maxtilt) ill_s],'double');
%% do fourier filtering for each illumination angle
schlieren_mask=rr(maxtilt)<max(size(maxtilt))/(4*zoom);
for a=0:ill_s(1)-1
    for b=0:ill_s(2)-1
        kx=double(ill_kx(a,b));
        ky=double(ill_ky(a,b));
        ill=exp(1i*2*pi*(kx.*xx(mma_phase)+ky.*yy(mma_phase)));
        u1=mma_phase.*ill;
        ftu1=ft(u1);
        vol_i2(:,:,a,b)=abs(ftu1).^2;
        u2=ftu1.*schlieren_mask;
        u3=ft(u2);
        i3=abs(u3).^2;
        vol_i3(:,:,a,b)=i3(:,:);
    end
end

%% display 3 fourier plane intensity images overlayed with aperture
schlieren_ring=(rr(maxtilt)>max(size(maxtilt))/(4*zoom))&...
    (rr(maxtilt)<max(size(maxtilt))*1.1/(4*zoom));
wide=horzcat(...
    1e7*normalize(squeeze(vol_i2(:,:,2,3))),...
    1e7*normalize(squeeze(vol_i2(:,:,6,6))),...
    1e7*normalize(squeeze(vol_i2(:,:,12,12))));
filt_im=overlay(wide,repmat(schlieren_ring,[3 1]));
writeim(filt_im,'01filt_im.jpg','JPEG');

%% display the corresponding intensity images of the mma
wide=horzcat(...
    normalize(squeeze(vol_i3(:,:,2,3))),...
    normalize(squeeze(vol_i3(:,:,6,6))),...
    normalize(squeeze(vol_i3(:,:,12,12))));
writeim(255*wide,'02mma_im.jpg','JPEG');

%% create a mosaic of all the images central image is illuminated
%% perpendicular
[w h a b]=size(vol_i3);
mosaic_i3=reshape(permute(reshape(vol_i3,[w h*b a]),[2 1 3]),[w*a h*b]);
mosaic_i3_small=resample(mosaic_i3,[1/zoom 1/zoom]);
mosaic_i2=reshape(permute(reshape(vol_i2,[w h*b a]),[2 1 3]),[w*a h*b]);
mosaic_i2_small=resample(mosaic_i2,[1/zoom 1/zoom]);

%% show (log of) intensity in the fourier plane with aperture
w=50; h=50;
ill_aps=reshape(permute( ...
    reshape(repmat(extract(schlieren_mask,[w h]),...
    [1 1 a b]),[w h*b a]),[2 1 3]),[w*a h*b]);
mosaic_i2_cut=reshape(permute(reshape(extract(vol_i2,...
    [w h]),[w h*b a]),[2 1 3]),[w*a h*b]);
mosaic_filter=overlay(1e5*normalize(mosaic_i2_cut),~ill_aps);
writeim(mosaic_filter,'03mosaic2.jpg','JPEG');

%% draw a ring into the mosaic. all images inside of the ring will be
%% incoherently summed to form the partial coherent image
ill_aperture=(rr(mosaic_i3_small,'freq')>.8/2) & ...
    (rr(mosaic_i3_small,'freq')<.84/2);
mosaic3_ring=overlay(255*mosaic_i3_small/max(mosaic_i3_small),...
    ill_aperture);
writeim(mosaic3_ring,'04mosaic3.jpg','JPEG');

%% do the incoherent summation
inco_accum=newim(zoom*s);
for a=0:ill_s(1)-1
    for b=0:ill_s(2)-1
        if ill_on(a,b)
            inco_accum(:,:)=inco_accum(:,:)+squeeze(vol_i3(:,:,a,b));
        end
    end
end
%%
writeim(255*normalize(inco_accum),'05incoherent','JPEG');
