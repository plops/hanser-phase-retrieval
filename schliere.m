% coherent image of mma after fourier filter

%%
zoom=5;
tm=(rr(s)<9.5); %.*mod(yy(s),2)
dir=2*(mod(yy(s),2)-.5); % mirror rows tilt in different directions
bigtilt=incimate((1-tm).*dir,zoom);

% for odd zoom, left edge of each mirror is -.25 and right edge +.25
shif=floor(zoom/2); 
if mod(zoom,2)==0
    shif=shif-1;
end
maxtilt=(mod(xx(bigtilt,'corner'),zoom)-shif)/(floor(zoom/2)*4);

mma_phase=exp(1i*2*pi*bigtilt.*maxtilt);

% kx in [-.5 .. .5]/zoom
kx=-.25/zoom; ky=.125/zoom;
ill=exp(1i*2*pi*(kx.*xx(mma_phase)+ky.*yy(mma_phase)));

u1=mma_phase.*ill;
u2=ft(u1).*(rr(u1)<max(size(u1))/(2*zoom));
u3=ft(u2);