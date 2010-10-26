% coherent image of mma after fourier filter

%%
zoom=5;
% tm=(rr(s)<4.5); %.*mod(yy(s),2) % hard edged circle
tm=exp(-((xx(s)-14).^2+yy(s).^2)/30); % soft gaussian
dir=2*(mod(yy(s),2)-.5); % mirror rows tilt in different directions
bigtilt=incimate((1-tm).*dir,zoom);

% for odd zoom, left edge of each mirror is -.25 and right edge +.25
shif=floor(zoom/2); 
if mod(zoom,2)==0
    shif=shif-1;
end
maxtilt=(mod(xx(bigtilt,'corner'),zoom)-shif)/(floor(zoom/2)*4);

mma_phase=exp(1i*2*pi*bigtilt.*maxtilt);

%%
ill_s=[13 13];
ill_kx=xx(ill_s,'freq')/(2*zoom);
ill_ky=yy(ill_s,'freq')/(2*zoom);
ill_on=rr(ill_s,'freq')<=.5;
%%

vol_i3=newim([size(maxtilt) ill_s],'double');
vol_i2=newim([size(maxtilt) ill_s],'double');
%%
% kx in [-.5 .. .5]/zoom
schlieren_mask=rr(u1)<max(size(u1))/(4*zoom);
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

%%
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
lmi2cut=log(mosaic_i2_cut);
overlay((lmi2cut-min(lmi2cut))*255/(max(lmi2cut)-min(lmi2cut)),~ill_aps)
%%
ill_aperture=(rr(mosaic_i3_small,'freq')>.8/2) & ...
    (rr(mosaic_i3_small,'freq')<.84/2);
overlay(255*mosaic_i3_small/max(mosaic_i3_small),ill_aperture)