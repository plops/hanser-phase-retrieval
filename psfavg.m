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
    n=neighbors(dt,i); % get all neighbors
    [rown,coln]=size(n);
    for j=1:coln  % distance to each of the neighbors
      if isnan(n(j))||n(j)>rowp % NaN when point with itself     
          nd(j)=1e12;
      else
          nd(j)=.04095*norm(p(i,:)-p(n(j),:));
      end
    end
    val(i)=min(nd); % minimal distance
    if val(i)==1e12
       val(i)=0.0;
    end
end
[sor,indsor]=sort(val);
% now we should select the last 5 points or so
np=4;
point=zeros(np,2);
for i=1:np
    point(i,:)=p(indsor(end-i),:);
end
dipshow(a),hold on
plot(point(:,1),point(:,2),'+')
