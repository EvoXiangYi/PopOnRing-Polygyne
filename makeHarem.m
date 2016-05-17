function [matrixWives] = makeHarem(femalesIdx, femalesFitness, malesIdx, malesFitness, meanNumWives)
	% Assume the number of wives a male have follows the poisson distribution with lambda = meanNumWives
	%generate cumulative distribution
	x = 0:2*meanNumWives;
	y = poisspdf(x,meanNumWives);
	cum_y = cumsum(y);
	cum_y(end) = 1; % a male can have 2 times meanNumWives at maximum.
	
	matrixWives=[];
	
	while and(length(femalesIdx)>0, length(malesIdx)>0)
		chosenMaleIdx = myID_parent(malesFitness);
		husbandIdx = malesIdx(chosenMaleIdx);
		actuallyNumWives=find(cum_y>rand,1)-1;% minus one because the number of wives start from 0.
		
		i=0; %count number of wives already taken
		while and(length(femalesIdx)>0, i<actuallyNumWives+1)
			i=i+1;
			chosenFemaleIdx = myID_parent(femalesFitness);
			wife_idx = femalesIdx(chosenFemaleIdx);
			wife_fitness = femalesFitness(chosenFemaleIdx);
			matrixWives=[matrixWives,[wife_idx; wife_fitness; husbandIdx]];
			
			femalesIdx(chosenFemaleIdx)=[];
			femalesFitness(chosenFemaleIdx)=[];
		end
		malesIdx(chosenMaleIdx)=[];
		malesFitness(chosenMaleIdx)=[];
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test the function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%femalesIdx = 1:2:51;
%femalesFitness = 30:-1:5;
%maleIdx=2:2:52;
%malesFitness = 30:-1:5;
%meanNumWives = 5;
%matrixWives = makeHarem(femalesIdx, femalesFitness, maleIdx, malesFitness, meanNumWives)