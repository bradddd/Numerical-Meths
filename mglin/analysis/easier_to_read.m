fid=fopen('soln.dat','r');   % open data file
fseek(fid,0,'eof');          % go to end of file
position = ftell(fid);       % determine size of file
frewind(fid);                % go back to the beginning
n = sqrt(position/8);        % assume file contains 2x2 matrix
f=fread(fid,[n,n],'double'); % of 8-byte doubles
fclose(fid);                 % close data file
surf(-f);                    % draw without contour lines
shading flat;                % disable black edges
colormap(hsv(128));          % use more and different colors
