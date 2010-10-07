%%
bg=readim('Snap.pgm');
a=readim('30.pgm');
a=a-bg;

%% find local maxima in regions that are bright enough to contain beads
a-a.*maxima(gaussf(a.*(a>200)))

%% find the coordinates of the maxima (even with subpixel accuracy)
p=findmaxima(gaussf(a.*(a>200)));
% overlay coordinates into image
dipshow(a),hold on
plot(p(:,1),p(:,2),'+')

%% construct delaunay triangulation (uses CGAL internally)
dt=DelaunayTri(p);
% overlay triangles on image
dipshow(a),hold on
triplot(dt)

%% find nearest neighbor for each of the points in p
for i=1:1
    n=neighbors(dt,i); % get all neighbors
    [p,q]=size(n);
    for j=1:q  % distance to each of the neighbors, ignore NaN
      if isnan(n(j))          
          nd(j)=1e12;
      else
          nd(j)=norm(p(i,:)-p(n(j),:));
      end
    end
    val(i)=min(nd); % minimal distance
end