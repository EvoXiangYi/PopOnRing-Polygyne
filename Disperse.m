function [HabitatsAfterDispersal, MigrationStatus]=Disperse(numHabitats, CurrentHabitats, ProbOffspringDispersal, DisperseDeathRate)
	
popSize=length(CurrentHabitats);
HabitatsAfterDispersal=zeros(1,popSize);
MigrationStatus = zeros(1, popSize);

for i=1:popSize
	currentHabitat=CurrentHabitats(i);
	willmove=rand<ProbOffspringDispersal(i);
	if willmove
		% Die from disperse, Habitat becomes 0
		if rand<DisperseDeathRate
			continue;
		end
		% Survived from dispersal
		if and(currentHabitat>1, currentHabitat<numHabitats)
			if rand<0.5
				currentHabitat=currentHabitat-1;
			else
				currentHabitat=currentHabitat+1;
			end
		elseif currentHabitat==1
			if rand<0.5
				currentHabitat=numHabitats;
			else
				currentHabitat=2;
			end
		else
			if rand<0.5
				currentHabitat=numHabitats-1;
			else
				currentHabitat=1;
			end
		end
		MigrationStatus(i)=1;
	end
	HabitatsAfterDispersal(i)=currentHabitat;
end

endfunction


%%%%%%%%%%%%%%%%% Test the Function %%%%%%%%%%%%%%%%%%%%%%%%%%%
%numHabitats=2;
%ProbOffspringDispersal=ones(1,10);
%DisperseDeathRate=0;
%Case of all disperse
%CurrentHabitats=[1,2,1,2,1,2,1,2,1,2]
%[NewHabitats,MigrantStatus]=Disperse(numHabitats, CurrentHabitats, ProbOffspringDispersal, DisperseDeathRate)
%Case of high death rates
%DisperseDeathRate=1;
%NewHabitats=Disperse(numHabitats, CurrentHabitats, ProbOffspringDispersal, DisperseDeathRate)
