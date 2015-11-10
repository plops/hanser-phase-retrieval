%% incoherent addition of images in focal plane
tic
s=[64 64];
tm=exp(-((xx(s)-14).^2+yy(s).^2)/30); 
tl=exp(-((xx(s)-2).^2+(yy(s)+4).^2)/100); % soft
tl=sqrt((xx(s)-12).^2+(yy(s)+20).^2)<10 | rr(s)<4;
angles=7;
%%
acc=newim([s (2*angles+1)^2+1],'dfloat');
ma=1/2;
i=0;
for kx=-ma:ma/angles:ma
    for ky=-ma:ma/angles:ma
        i=i+1;
        acc(:,:,i)=squeeze(acc(:,:,i-1))+...
            abs(memi(kx,ky,tm,tl)).^2;
    end
end
toc


