function [fitnessVector]=FitnessLite(envTraitVector, habitatCenter, fitnessIOS)
	popSize=length(envTraitVector);
	distances = abs(envTraitVector - habitatCenter);
	fitnessVector = exp(- fitnessIOS * distances .^ 2);
end


%%%%%%%%%%%%%%%%%%%%%% Test the Function %%%%%%%%%%%%%%%%%%%%%%%%%
%envTraitVector=0.1:0.1:0.9;
%habitatCenter=0.5;
%fitnessIOS=3;
%fitnessVector=FitnessLite(envTraitVector, habitatCenter, fitnessIOS)