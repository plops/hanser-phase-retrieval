%% incoherent addition of images in focal plane
tic
s=[64 64];
angles=7;
acc=newim([s (2*angles+1)^2+1],'dfloat');
ma=1/2;
i=0;
for kx=-ma:ma/angles:ma
    for ky=-ma:ma/angles:ma
        i=i+1;
        acc(:,:,i)=squeeze(acc(:,:,i-1))+...
            abs(memi(kx,ky,rr(s)<8.3,rr(s)<30)).^2;
    end
end
toc


