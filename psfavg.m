% select regions in a bead stack that contain single beads
%%
data='/home/martin/d1004/normal_0/pgm/'
bg=readim([data 'Snap.pgm']);
a=readim([data '30.pgm']);
a=a-bg;

%% find the coordinates of the maxima (even with subpixel accuracy)
p=findmaxima(gaussf(a.*(a>400)));
% overlay coordinates into image
dipshow(a),hold on
plot(p(:,1),p(:,2),'+')

%% construct delaunay triangulation (uses CGAL internally)
dt=DelaunayTri(p);
% overlay triangles on image
dipshow(a),hold on
triplot(dt)

%% find nearest neighbor for each of the points in p
[rowp,colp]=size(p);
for i=1:rowp
    n=cell2mat(dt.vertexAttachments(i)); % get all neighboring triangles
    [rown,coln]=size(n);
    val(i)=1e12;
    for j=1:coln  % check each neighbor triangle
        for k=1:3 % check each point on the triangle %.04095
            if dt(n(j),k)==i
            else
                val(i)=min(val(i),norm(p(i,:)-p(dt(n(j),k),:)));
            end
        end
    end
end
[sor,indsor]=sort(val);
% now we should select the last 5 points that are not too close to border
border=60;
np=12;
point=zeros(np,2);
i=0;
k=1;
while k<=np
    j=indsor(end-i);
    % don't consider points that are to near to the border
    if border<p(j,1) && p(j,1)<1391-border && ...
       border<p(j,2) && p(j,2)<1039-border
        point(k,:)=p(j,:);
        k=k+1;
    end
    i=i+1;
end
dipshow(a),hold on
plot(point(:,1),point(:,2),'x','MarkerSize',22,'LineWidth',4,'MarkerEdgeColor','g')
triplot(dt)
for k=1:rowp
    text(p(k,1),p(k,2),sprintf('%3.0f',val(k)),'Color','r','FontSize',18)
end
%% read in the volume
% the first image is rubbish (z-stage moved during exposure)
% z-step is 100nm, the focus is approximately at i=32
% in x and y .04095 um/pixel 
% (6.45/(2.5*63) pixel pitch 6.45um, optovar 2.5, objmag 63)
nimg=61;
[sx,sy]=size(bg);
stack=newim(sx,sy,nimg,datatype(bg));
for i=0:nimg-1
    i
    stack(:,:,i)=readim([data sprintf('%02d.pgm',i)])-bg;
end

%% copy each lonely point into an array of psfs
dimsh=[120 120 60];
psfs=newim([dimsh np]);
for k=1:np
    psfs(:,:,:,k-1)=extract(stack,dimsh,[floor(point(k,:)) 32]);
end

%%
corrected=psfs;
a=DampEdge(squeeze(psfs(:,:,:,0)),0.8,3,0);
sv=zeros(np,3);
for k=0:np-1
    b=DampEdge(squeeze(psfs(:,:,:,k)),0.8,3,0);
    sv(1+k,:)=findmaxima(abs(ift(ft(a).*conj(ft(b)))),'gaussian')-floor(dimsh/2)
    %sv(1+k,:)=findshift(a,b,'iter')
    corrected(:,:,:,k)=shift(squeeze(psfs(:,:,:,k)),sv(1+k));
end

% findshift
%     0.0000    0.0000   -0.0000
%    -0.2149   -0.3264    0.6406
%    -0.0258   -0.1172    0.6919
%    -0.7813   -0.2322    0.7846
%    -0.7335   -0.1299    0.0075
%    -0.2322   -0.6418    0.2194
% my cross correlation
%     0         0         0
%    -0.2955   -0.4740    0.5145
%    -0.0192   -0.1492    0.5799
%    -0.6683   -0.3284    0.6925
%    -0.6683   -0.1677    0.0070
%    -0.2246   -0.6782    0.1989
%%
psf=squeeze(mean(corrected,[],4))
otf=ft(psf)
mean(std(corrected,[],4)) % FIXME maybe individual scaling is necessary
otf2d=mean(otf,[],3)
%otf_alternate=mean(dip_fouriertransform(psfs,'forward',[1 1 1 0]),[],4)
%%
for k=[-25,-20,-15,0,15,20,25]
    writeim(squeeze(psf(:,:,29+k)),...
        sprintf('/home/martin/1007/avgpsf/psf%+03d',k),'ICSv2');
end
%% rotational integration
nr=100;
otf1d=newim(nr);
for r=0:nr-1
    rad=r/(nr+1)*max(dimsh(1:2))/2;
    nphi=r+1;
    sum=0;
    for phii=0:nphi
        phi=2*pi/nphi*phii;
        sum=sum+...
        interp2(double(abs(otf2d)),...
            floor(dimsh(1)/2)+rad*cos(phi),...
            floor(dimsh(2)/2)+rad*sin(phi),'cubic');
    end
    otf1d(r)=sum/nphi;
end
otf1d