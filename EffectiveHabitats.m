function [effectiveHabitats]=EffectiveHabitats(numHabitats,habitatVector,sexVector)
	effectiveHabitats=[];
	for i=1:numHabitats
		if and(sum(and(habitatVector==i,sexVector==0))!= 0, sum(and(habitatVector==i,sexVector==1))!= 0)
			effectiveHabitats=[effectiveHabitats,i];
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test the Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%numHabitats=4;
%habitatVector=[0,0,1,2,3,4,0,0,1,2,3,4,0,0,1,2,3,4,0,0,1,2,3,4]
%sexVector=[0,0,1,0,0,0,0,0,1,0,0,0,1,1,1,1,1,0,1,1,1,1,1,0]
%effectiveHabitats=EffectiveHabitats(numHabitats,habitatVector,sexVector)