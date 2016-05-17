function [carryingCapacitySurvivalVector]=CarryingCapacityRegulation(juvenileHabitatsVector, numHabitats, carryingCapacity)
	popSize=length(juvenileHabitatsVector);
	carryingCapacitySurvivalVector=zeros(1,popSize);
	for i=1:numHabitats
		Juv_on_Habitat_i=find(juvenileHabitatsVector==i);
		num_Juv_on_Habitat_i=length(Juv_on_Habitat_i);
		if or(num_Juv_on_Habitat_i < carryingCapacity, num_Juv_on_Habitat_i == carryingCapacity)
			carryingCapacitySurvivalVector(Juv_on_Habitat_i)=1;
		else
			carryingCapacitySurvivalVector(Juv_on_Habitat_i(randperm(num_Juv_on_Habitat_i)(1:carryingCapacity)))=1;
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test the Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%juvenileHabitatsVector=[1,2,3,1,2,3,1,2,3,1,1,1,1,1,1,1]
%numHabitats=3;
%carryingCapacity=3;
%surVector=CarryingCapacityRegulation(juvenileHabitatsVector, numHabitats, carryingCapacity)