matsize = 10; 
p=matsize;
q=matsize;
r=matsize;

fp = fopen('C:\Octave\output_CN_10x10x10.bin');


B = fread(fp,Inf,'double');

fclose(fp);

size(B)

nummatrices = size(B)/(matsize^3)

Bans = reshape(B, r, q, p, nummatrices);

for t=1:nummatrices
	surfc(Bans(:,:,ceil((q)/2),t));
	set(gca,'ZLim',[0,5]);
	axis square;
	%pause(.005);
	pause(5);
end


%surf(Bans(:,:,:))


%Bans = reshape(B, r, q, p);
%surfc(Bans(:,:,10))