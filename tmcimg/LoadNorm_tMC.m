% function [a, bgs, fluence, pOpt] = LoadNorm_tMC( filenm )


function [a, bgs, fluence, pOpt] = LoadNorm_tMC( filenm )

fid=fopen(sprintf('%s.inp',filenm),'r');
NPh = fscanf(fid,'%d', 1);
foo = fscanf(fid,'%d', 1);
xi = fscanf(fid,'%d',1 );
yi = fscanf(fid,'%d',1 );
zi = fscanf(fid,'%d',1 );
cxi = fscanf(fid,'%d',1 );
cyi = fscanf(fid,'%d',1 );
czi = fscanf(fid,'%d',1 );
minT = fscanf(fid,'%f',1 );
maxT = fscanf(fid,'%f',1 );
stepT = fscanf(fid,'%f',1 );   
binFile = fscanf(fid,'%s',1 );
xstep = fscanf(fid,'%f', 1);
nxstep = fscanf(fid,'%d', 1);
xmin = fscanf(fid,'%d',1 );
xmax = fscanf(fid,'%d',1 );
ystep = fscanf(fid,'%f', 1);
nystep = fscanf(fid,'%d', 1);
ymin = fscanf(fid,'%d',1 );
ymax = fscanf(fid,'%d',1 );
zstep = fscanf(fid,'%f', 1);
nzstep = fscanf(fid,'%d', 1);
zmin = fscanf(fid,'%d',1 );
zmax = fscanf(fid,'%d',1 );
nTissue = fscanf(fid,'%d',1 );
for idx=1:nTissue
  tis(idx,:) = fscanf(fid,'%f',4 )';
end
nDet = fscanf( fid, '%d', 2);
for idx=1:nDet(1)
  pOpt(:,idx) = fscanf( fid, '%d', 3 );
end

fclose(fid);


%xmin = xmin + 1;
%ymin = ymin + 1;
%zmin = zmin + 1;

%xmax = xmax + 1;
%ymax = ymax + 1;
%zmax = zmax + 1;

nx = xmax-xmin+1;
ny = ymax-ymin+1;
nz = zmax-zmin+1;

fid=fopen(sprintf('%s.2pt',filenm),'r');
a = fread(fid,'float64');
fclose(fid);
nt = length(a) / (nx*ny*nz);
a = reshape(a,[nx,ny,nz,nt]);
a = a/NPh;

fid=fopen(binFile,'rb');
bg = fread(fid,nxstep*nystep*nzstep,'int8');
fclose(fid);
bg = reshape(bg,[nxstep,nystep,nzstep]);
bgs = bg(xmin:xmax, ymin:ymax, zmin:zmax );

b=sum(a,4);
for idx=1:nTissue
  list = find( bgs==idx );
  b(list) = b(list) * tis(idx,3);
end

list = find(a<0);
am = sum(a(list));

list = find(b>0);
ap = sum(b(list));

list = find(a>0);

ak = (1+am) / ap;
a(list) = ak * a(list);

if nt>1
  a = a / stepT;
end

if nDet(1)>0
  pOpt(1,:) = pOpt(1,:) - xmin + 1;
  pOpt(2,:) = pOpt(2,:) - ymin + 1;
  pOpt(3,:) = pOpt(3,:) - zmin + 1;
end

fluence = zeros(nDet(1),1);
for idx=1:nDet(1)
  fluence(idx) = a(pOpt(1,idx), pOpt(2,idx), pOpt(3,idx));
end

  