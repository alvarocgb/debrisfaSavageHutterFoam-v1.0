forAll(boundaryNames, boundaryI)
{
	const word boundaryName = boundaryNames[boundaryI];

	dictionary subHidrogram = Hidrograms.subDict(boundaryName);

	List<scalar> subwide  = subHidrogram.lookup("wide");
	List<scalar> subslope = subHidrogram.lookup("slope");
	List<vector> subdirections = subHidrogram.lookup("directions");
	List<scalar> subheight = subHidrogram.lookup("heights");	
	List<scalar> subtimes  = subHidrogram.lookup("times");
	List<scalar> subflows  = subHidrogram.lookup("flows");
	List<scalar> subCvs = subHidrogram.lookup("Cvs");	
	scalar alpha = subHidrogram.lookupOrDefault<scalar>("alpha",0);

	label t = runTime.timeOutputValue();

	float Qt  = 0;
	float Cvt = 0;
	label i = 0;
	label j = subtimes.size()-1;

	if (t<subtimes[i])
        {
            FatalIOErrorIn("createFlowCondition.H", releaseFlow)
                << "Time " << runTime.timeName() << "s is not a"
                << " valid time for the introduced hidrograph." << exit(FatalIOError);
        }
	else if (t<= subtimes[j])
	{
		while(j-i>1)
		{
			int n = i+(j-i)/2;
			if (subtimes[n]==t)
			{
				i = n;
				j = n;
			}
			else if (subtimes[n]>t)
			{
				j = n;
			}
			else
			{
				i = n;
			}
        }
		Qt = subflows[i]+(subflows[j]-subflows[i])/(subtimes[j]-subtimes[i]+SMALL)*(t-subtimes[i]);
		Cvt = subCvs[i]+(subCvs[j]-subCvs[i])/(subtimes[j]-subtimes[i]+SMALL)*(t-subtimes[i]);
	}
	else
	{
		Qt = subflows[j];
		Cvt = subCvs[j];
	}

	float rhot = (1-Cvt)*rho_w.value()+Cvt*rho_s.value();
	float Cwt  = 1-Cvt;

	const scalar & ut = friction->ut(Qt*alpha, subwide, subslope, rhot, boundaryName);	
        const List<scalar> & ht = friction->ht(Qt*alpha, subwide, subslope, subheight, rhot, boundaryName);

	label length(faBoundaries[faboundaryNames.find(boundaryName)].size());
	label index(faboundaryNames.find(boundaryName));

	areaScalarField::Boundary boundaryh(h.boundaryFieldRef());
	Field<scalar> field_ht(ht);
	boundaryh[index] == field_ht;
	h.boundaryFieldRef() == boundaryh;

	areaVectorField::Boundary boundaryUs(Us.boundaryFieldRef());
	Field<vector> field_ut(subdirections*ut);
	boundaryUs[index] == field_ut;
	Us.boundaryFieldRef() == boundaryUs;

	areaScalarField::Boundary boundaryrho(rho.boundaryFieldRef());
	Field<scalar> field_rhot(length, rhot);
	boundaryrho[index] == field_rhot;
	rho.boundaryFieldRef() == boundaryrho;

	areaScalarField::Boundary boundaryrhoh(rhoh.boundaryFieldRef());
	Field<scalar> field_rhoht(rhot*ht);
	boundaryrhoh[index] == field_rhoht;
	rhoh.boundaryFieldRef() == boundaryrhoh;

	areaScalarField::Boundary boundaryCw(Cw.boundaryFieldRef());
	Field<scalar> field_Cwt(length, Cwt);
	boundaryCw[index] == field_Cwt;
	Cw.boundaryFieldRef() == boundaryCw;

	areaScalarField::Boundary boundaryCwh(Cwh.boundaryFieldRef());
	Field<scalar> field_Cwht(Cwt*ht);
	boundaryCwh[index] == field_Cwht;
	Cwh.boundaryFieldRef() == boundaryCwh;
};
