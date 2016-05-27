function []= DispersalSimulation_Polygyne(InMeanNumWives, InTimeSteps, InNumGrid, InPopSize, InNumHabitats, InpSpatial, InpTemporal, InfitnessIOS, InMeanDispersalTrait_init, InSigmaDispersalTrait_init, InSigmaTraitCoordinates, InSigmaDispersal, InDispersalDeathRate, InNumRepeats)
	
	%%%%%%%%%%%%%%%%%%%%% start Parameters Initialization  %%%%%%%%%%%%%%%%%%%%%%%%%
	meanNumWives = InMeanNumWives;
	numRepeats = InNumRepeats;
	timeSteps = InTimeSteps;
	numGrid = InNumGrid;
	popSize = InPopSize;
	numHabitats = InNumHabitats;
	carryingCapacity = 2 * popSize/numHabitats; %carrying capacity of each of the habitats
	pSpatial=InpSpatial; % Spatial autocorrelation parameter
	pTemporal=InpTemporal; % Temporal autocorrelation parameter
	fitnessIOS=InfitnessIOS;
	init_Mean_Dispersal = InMeanDispersalTrait_init;
	init_Prob_dispersal=1/(1+exp(-init_Mean_Dispersal)); %probability of dispersal for both males and females
	init_Sigma_Dispersal = InSigmaDispersalTrait_init;
	sigmaTraitCoordinates=InSigmaTraitCoordinates;
	sigmaDispersal=InSigmaDispersal;
	dispersalDeathRate=InDispersalDeathRate;
	%%%%%%%%%%%%%%%%%%%%%  end Parameters Initialization  %%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%%%%%%%%%%%%%%%%%%%%  start Habitats Initialization  %%%%%%%%%%%%%%%%%%%%
	%To always start from the same environment
	randomStete=csvread("randStateFile.csv");
	rand("state",randomStete);
	inputEnvironment=rand(1,numGrid); 	%In Luke and Hanna's paper: inputEnvironment=unidrnd(2,1,numGrid)-1;
	envrionmentVector=UpdateEnvironment(inputEnvironment,numGrid,pSpatial,pTemporal,500);
	habitatLocations=SetHabitats(numGrid, numHabitats);
	habitatCenters=envrionmentVector(habitatLocations);
	%%%%%%%%%%%%%%%%%%%%%  end Habitats Initialization  %%%%%%%%%%%%%%%%%%%%%%
	
	
	%%%%%%%%%%%%%%%%%%%%%  start GenomeStructure  %%%%%%%%%%%%%%%%%%%%%%
	envTraitLocus=1:2; %adaptation traits from both mom and dad
	dispLociF=3:4;
	dispLociS=5:6;
	dispLociM=7:8;
	sex=9;
	habitat=10;
	migrateStatus = 11;
	numLoci=migrateStatus;
	%%%%%%%%%%%%%%%%%%%%%  end GenomeStructure  %%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%%%%%%%%%%%%%%%%%%%%  start PopulationStructure  %%%%%%%%%%%%%%%%%%%%
	population=zeros(numLoci, popSize); %Initialize the "diploid genome"
	population(envTraitLocus,:) = rand(envTraitLocus(end), popSize); %Initialize environmatal coordinates
	population(dispLociF(1):dispLociM(end),:) = normrnd(init_Mean_Dispersal,init_Sigma_Dispersal,[dispLociM(end)-envTraitLocus(end),popSize]);
	population(sex,:) = unidrnd(2,1,popSize)-1; %female:0, males:1
	population(habitat,:) = unidrnd(numHabitats,1,popSize);
	population(migrateStatus,:) = 0; %Starting from no migrants
	%%%%%%%%%%%%%%%%%%%%%  end PopulationStructure  %%%%%%%%%%%%%%%%%%%%%%
	
	
	%%%%%%%%%%%%%%%%%%%%%  start Output Data Structure  %%%%%%%%%%%%%%%%%%%%
	% 1->Total number of individuals
	% 2->Total numbers of Habitats that have both males and females
	% 3->Mean fitness of native females
	% 4->Mean fitness of migrant females
	% 5->Mean fitness of native males
	% 6->Mean fitness of migrant males
	% 7->Mean dispersal probability for females
	% 8->Mean dispersal probability for males
	SummaryData=NaN(8,timeSteps);
	%%%%%%%%%%%%%%%%%%%%%  end Output Data Structure  %%%%%%%%%%%%%%%%%%%% 
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%% start Population Update    %%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for t=1:timeSteps
		t
		
		%%%%%%%%%%%%%%%%%%%%%  start Environment Selection  %%%%%%%%%%%%%%%%%%
	
		%%%%%%%%%%%%%%%%%%%%%  end Environment Selection  %%%%%%%%%%%%%%%%%%%%
	
	
		%%%%%%%%%%%%%%%%%%%%%  start Carrying Capacity regulation  %%%%%%%%%%%%%%%%%%
		surviveVectorAfterRegulation = CarryingCapacityRegulation(population(habitat,:), numHabitats, carryingCapacity);
		population(habitat,surviveVectorAfterRegulation==0)=0;
		%%%%%%%%%%%%%%%%%%%%%  end Carrying Capacity regulation  %%%%%%%%%%%%%%%%%%%%
	
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  start population census  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%% 1 ->Total number of individuals : totalNumIndividuals
		totalNumIndividuals=sum(population(habitat,:)!=0);
		if totalNumIndividuals < 1
			break; %exit simulation upon extinction.
		end
		
		%% 2 ->Total numbers of Habitats that have both males and females : totalEffectiveHabitats
		effectiveHabitats=EffectiveHabitats(numHabitats,population(habitat,:),population(sex,:));
		totalEffectiveHabitats=length(effectiveHabitats);
		
		%% 3,4,5,6 ->Mean fitness for native females and total females, mean fitness for native males and total males : meanFemaleFitnessNative, meanFemaleFitness, meanMaleFitnessNative, meanMaleFitness
		femalesOnHabitatsCell=cell(1,totalEffectiveHabitats);
		femalesFitnessOnHabitatsCell=cell(1,totalEffectiveHabitats);
		malesOnHabitatsCell=cell(1,totalEffectiveHabitats);
		malesFitnessOnHabitatsCell=cell(1,totalEffectiveHabitats);
		
		NativeFemales = [];
		NativeFemalesFitness = [];
		NativeMales = [];
		NativeMalesFitness = [];
		
		for i=1:totalEffectiveHabitats
			localHabitat=effectiveHabitats(i);
			%% Distinguish between native and migrant females
			% find indices
			idx_femalesN=find(and(population(habitat,:)==localHabitat, population(sex,:)==0, population(migrateStatus,:)==0));
			NativeFemales = [NativeFemales, idx_femalesN];
			
			idx_females=find(and(population(habitat,:)==localHabitat, population(sex,:)==0));
			femalesOnHabitatsCell{i} = idx_females;
			
			% find fitness values
			femaleFitnessN = FitnessOrder1(mean(population(envTraitLocus,idx_femalesN)),habitatCenters(localHabitat),fitnessIOS);
			NativeFemalesFitness = [NativeFemalesFitness, femaleFitnessN];
			
			femaleFitness = FitnessOrder1(mean(population(envTraitLocus,idx_females)),habitatCenters(localHabitat),fitnessIOS);
			femalesFitnessOnHabitatsCell{i} = femaleFitness;
			
			%% Distinguish between native and migrant males
			% find indices
			idx_malesN = find(and(population(habitat,:)==localHabitat, population(sex,:)==1, population(migrateStatus,:)==0));
			NativeMales = [NativeMales, idx_malesN];
			
			idx_males = find(and(population(habitat,:)==localHabitat, population(sex,:)==1));
			malesOnHabitatsCell{i} = idx_males;
			
			% find fitness values
			maleFitnessN = FitnessOrder1(mean(population(envTraitLocus,idx_malesN)),habitatCenters(localHabitat),fitnessIOS);
			NativeMalesFitness = [NativeMalesFitness, maleFitnessN];
			
			maleFitness = FitnessOrder1(mean(population(envTraitLocus,idx_males)),habitatCenters(localHabitat),fitnessIOS);
			malesFitnessOnHabitatsCell{i} = maleFitness;
		end
		
		meanFemaleFitnessNative = mean(NativeFemalesFitness);
		meanFemaleFitness = mean(cell2mat(femalesFitnessOnHabitatsCell));
		meanMaleFitnessNative = mean(NativeMalesFitness);
		meanMaleFitness = mean(cell2mat(malesFitnessOnHabitatsCell));
		
		%% 7,8 ->Mean dispersal probability for females and males : meanProbFemaleDispersal, meanProbMaleDispersal
		idx_females = cell2mat(femalesOnHabitatsCell);
		femaleDispersalTraitVector = mean(population(dispLociF(1):dispLociS(end),:))(idx_females);
		femaleDispersalProbability = 1./(1+exp(-femaleDispersalTraitVector));
		meanProbFemaleDispersal = mean(femaleDispersalProbability);
		
		idx_males = cell2mat(malesOnHabitatsCell);
		maleDispersalTraitVector = mean(population(dispLociS(1):dispLociM(end),:))(idx_males);
		maleDispersalProbability = 1./(1+exp(-maleDispersalTraitVector));
		meanProbMaleDispersal = mean(maleDispersalProbability);
		
		%% Record in SummaryData
		SummaryData(:,t)=[totalNumIndividuals; totalEffectiveHabitats; meanFemaleFitnessNative; meanFemaleFitness; meanMaleFitness; meanMaleFitnessNative; meanProbFemaleDispersal; meanProbMaleDispersal];
		%%%%%%%%%%%%%%%%%%%%%  end population census  %%%%%%%%%%%%%%%%%%%%%%%
		
		
		%%%%%%%%%%%%%%%%%%%%% getting married %%%%%%%%%%%%%%%%%%%%%
		MatrixWife=[];
		% row 1 -> indices of wives
		% row 2 -> fitness of wives
		% row 3 -> indices of husbands
		
		for i=1:totalEffectiveHabitats
			newWives=makeHarem(femalesOnHabitatsCell{i}, femalesFitnessOnHabitatsCell{i}, malesOnHabitatsCell{i}, malesFitnessOnHabitatsCell{i}, meanNumWives);
			MatrixWife=[MatrixWife, newWives];			
		end
		
		
		%%%%%%%%%%%%%%%%%%%%%%%% picking reproducing pairs %%%%%%%%%%%%%%%%%%%
		ParantMatrix = storkPicks(popSize, MatrixWife(1,:), MatrixWife(2,:), MatrixWife(3,:));
		
		%%%%%%%%%%%%%%%%%%%%%%%% making babies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		babyMatrix = makeBabies(population, ParantMatrix(1,:), ParantMatrix(2,:));
		
		population = babyMatrix;
		
		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dispersing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		traitsOffspringDispersal = mean(population(dispLociF(1):dispLociS(end),:)).*(population(sex,:)==0)+mean(population(dispLociS(1):dispLociM(end),:)).*(population(sex,:)==1);
		probOffspringDispersal = 1./(1+exp(-traitsOffspringDispersal));
		[population(habitat,:), population(migrateStatus,:)]= Disperse(numHabitats, population(habitat,:), probOffspringDispersal, dispersalDeathRate);
		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% updating environment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		envrionmentVector=UpdateEnvironment(envrionmentVector,numGrid,pSpatial,pTemporal,1);
		habitatCenters=envrionmentVector(habitatLocations);
	
	end
	
	%%%%%%%%%%%%%%%%%%%% Write SummaryData to Output %%%%%%%%%%%%%%%%%%%%%%%%%%%
	csvwrite([strcat("PopOnRing-Polygyne-InitDispProb-", num2str(init_Prob_dispersal), "-numGrid-", num2str(numGrid), "-numHabitats-", num2str(numHabitats),"-Rep-", num2str(numRepeats), ".csv")], SummaryData);
end