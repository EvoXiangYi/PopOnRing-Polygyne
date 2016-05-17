function [HabitatLocations]=SetHabitats(numGrid, numHabitats)
	x=1+round(linspace(0,numGrid,numHabitats+1));
	HabitatLocations=x(1:end-1);
endfunction

%%%%%%%%%%%%Test of function%%%%%%%%%%%%%%%%%%%%%%%
%numGrid=200;
%numHabitats=4;
%HabitatLocations=SetHabitats(numGrid, numHabitats)