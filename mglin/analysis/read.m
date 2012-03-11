n=1025;
fid=fopen('soln.dat','r');
f=fread(fid,[n,n],'double');
fclose(fid);
surfc(-f);
%contour(-f,10);