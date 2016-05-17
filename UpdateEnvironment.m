function [EnvrionmentVector]=UpdateEnvironment(OldEnvironment,numGrid,pSpatial,pTemporal,numUpdates)
	EnvrionmentVector=OldEnvironment;
	for tUpdates=1:numUpdates
		PreviousEnvionment=EnvrionmentVector;
		for i=1:numGrid
			%find direct neighbours and over border values
			ni=[i-1, i, i+1, numGrid+i-1, i+1-numGrid];
			%remove over border values
			f=find(ni>numGrid | ni<1); if ~isempty(f) ni(f)=[]; end;
			%mean values of the neighbouring cells
			Mean_correlated=mean(EnvrionmentVector(ni));
			%Update Habitat trait values
			EnvrionmentVector(i)=(pSpatial * Mean_correlated + (1 - pSpatial) * rand(1));
		end
		EnvrionmentVector = pTemporal * PreviousEnvionment + (1-pTemporal) * EnvrionmentVector;
	end
endfunction

%%%%%%%%%%%%Test of function%%%%%%%%%%%%%%%%%%%%%%%
%numGrid=20;
%pSpatial=0.9;
%pTemporal=0.9;
%InputEnvironment=unidrnd(2,1,numGrid)-1;
%OldEnvironment=UpdateEnvironment(InputEnvironment,numGrid,pSpatial,pTemporal,1000)
%NewEnvironment=UpdateEnvironment(OldEnvironment,numGrid,pSpatial,pTemporal,1)
