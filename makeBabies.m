%%%%%%%%%%%%%%%%%%%%%  IMPORTANT!!! CHECK Genome Structure  %%%%%%%%%%%%%%%%%%%%%%
%envTraitLocus=1:2; %adaptation traits from both mom and dad
%dispLociF=3:4;
%dispLociS=5:6;
%dispLociM=7:8;
%sex=9;
%habitat=10;
%migrateStatus = 11;
%numLoci=migrateStatus;


function [babyMatrix] = makeBabies(populationMatrix, MomIDs, DadIDs)
	popSize = size(populationMatrix)(2);
	babyMatrix = zeros(size(populationMatrix));
	MomMatrix = populationMatrix(:,MomIDs);
	DadMatrix = populationMatrix(:,DadIDs);
	
	momAlleles = [];
	for i=1:4 %4 is the number of diploid loci
		momAlleles = [momAlleles; pickOneAlleleAtLocus(popSize)];
	end
	
	dadAlleles = 1 - momAlleles;
	babyMatrix(1:8,:) = MomMatrix(1:8,:) .* momAlleles + DadMatrix(1:8,:) .* dadAlleles;
	babyMatrix(9,:) = unidrnd(2,1,popSize)-1;
	babyMatrix(10,:) = MomMatrix(10,:);
end

function[myAllele] = pickOneAlleleAtLocus(popSize)
	myRow = unidrnd(2,1,popSize)-1;
	myAllele = [myRow; 1-myRow];
end

%%%%%%%%%%%%%%%%%%% Test the function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%populationMatrix = [1,2,3,4,5].*ones(11,5);
%MomIDs = [1,2,2,3,3]
%DadIDs = [4,4,4,5,5]
%babyMatrix = makeBabies(populationMatrix, MomIDs, DadIDs)