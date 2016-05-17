function [id_parent] = myID_parent(fitnessVector)
	CDFVector = cumsum(fitnessVector);
	normCDFVector = CDFVector/CDFVector(end);
	x = rand;
	id_parent = find(normCDFVector>x,1);
endfunction