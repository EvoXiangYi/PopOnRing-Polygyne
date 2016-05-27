function [fitnessVector]=FitnessOrder1(envTraitVector, habitatCenter, fitnessIOS)
	%assume the maximum number off offspring a female can produce is 10;
	popSize=length(envTraitVector);
	distances = abs(envTraitVector - habitatCenter);
	fitnessVector = exp(- fitnessIOS * distances);
end


%%%%%%%%%%%%%%%%%%%%%% Test the Function %%%%%%%%%%%%%%%%%%%%%%%%%
%envTraitVector=0.1:0.1:0.9;
%habitatCenter=0.5;
%fitnessIOS=3;
%fitnessVector=FitnessLite(envTraitVector, habitatCenter, fitnessIOS)