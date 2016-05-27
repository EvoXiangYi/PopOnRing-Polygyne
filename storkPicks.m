function [parentMatrix] = storkPicks(popSize, femalesIdx, femalesFitness, malesIdx)
	parentPairs = [];
	for i=1:length(femalesIdx)
		newParentPairs = [femalesIdx(i); malesIdx(i)] .* ones(1,1+poissrnd(10*femalesFitness(i)));
		parentPairs = [parentPairs, newParentPairs];
	end
	storkPickedPairs = randperm(size(parentPairs)(2),popSize);
	parentMatrix = parentPairs(:,storkPickedPairs);
end


%%%%%%%%%%%%%%%%%%%%%% Test my function %%%%%%%%%%%%%%%%%%%%%%%%%%%
%popSize = 5;
%femalesIdx = [1,3,5,7,9];
%femalesFitness = [1,3,5,7,9];
%malesIdx = [2,4,6,8,10];
%PickedMatrix = storkPicks(popSize, femalesIdx, femalesFitness, malesIdx)