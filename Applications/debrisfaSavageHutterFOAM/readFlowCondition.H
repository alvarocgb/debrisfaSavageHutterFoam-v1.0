Info << nl << "Reading flow condition " << endl;
IOdictionary releaseFlow
(
    IOobject
    (
        "releaseFlow",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

IOdictionary gHidrograms
(
    IOobject
    (
        "Hidrograms",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const PtrList<entry> boundaries
(
    releaseFlow.lookup("boundaries")
);

wordList boundaryNames;
boundaryNames.setSize(boundaries.size());

faBoundaryMesh faBoundaries
(
    IOobject
    (
        "faBoundary",
        runTime.constant(),
		"faMesh",
		mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    aMesh
);

wordList faboundaryNames = faBoundaries.names();
wordList faboundaryTypes = faBoundaries.types();

dictionary Hidrograms("hidrograms");

forAll(boundaryNames, boundaryI)
{
	const entry& boundaryInfo = boundaries[boundaryI];
	
	if (!boundaryInfo.isDict())
	{
		FatalIOErrorIn("readFlowCondition.H", releaseFlow)
			<< "Entry " << boundaryInfo << " in boundaries section is not a"
			<< " valid dictionary." << exit(FatalIOError);
	}

	boundaryNames[boundaryI] = boundaryInfo.keyword();

	const word boundaryName = boundaryNames[boundaryI];

	if (faboundaryNames.find(boundaryName) == -1)
	{
		FatalIOErrorIn("readFlowCondition.H", releaseFlow)
			<< "Entry " << boundaryName << " in boundaries section is not a"
			<< " valid boundary in faBoundary." << exit(FatalIOError);
	}

    dictionary boundaryDict = boundaryInfo.dict();
	
	word sense;
	word direction;

	boundaryDict.lookup("sense") >> sense;
	boundaryDict.lookup("direction") >> direction;

    	keyType name = boundaryNames[boundaryI];

	dictionary subHidrogram(name);

	vector vector_dir(0,0,0);

	if (direction == "x")
	{
		vector_dir.x() = 1;
	}
	else if (direction == "y")
	{
		vector_dir.y() = 1;
	}
	else if (direction == "z")
	{
		vector_dir.z() = 1;
	}
	else
	{
		FatalIOErrorIn("readFlowCondition.H", releaseFlow)
			 << direction << " is not a valid Direction." << exit(FatalIOError);
	}

	scalar sign = 0;

	if (sense == "positive")
	{
		sign = 1;
	}
	else if (sense == "negative")
	{
		sign = -1;
	}
	else
	{
		FatalIOErrorIn("readFlowCondition.H", releaseFlow)
			<< sense << " is not a valid Sense." << exit(FatalIOError);
	}

	label index = faBoundaries.findIndex(boundaryName);

	const faPatch faPatchDict(faBoundaries[index],faBoundaries);	
	const vectorField edgevectorlength = faPatchDict.edgeLengths();
	const scalarField edgelength = faPatchDict.magEdgeLengths();
	const vectorField edgecentres = faPatchDict.edgeCentres();

	List<scalar> wide;
	List<scalar> slope;
	List<vector> directions;
	List<scalar> heights;
	scalar sumwide = 0;

	forAll(edgelength, edgeI)
	{
		wide.append(edgelength[edgeI]);
		sumwide += edgelength[edgeI];
		vector unitary_vector = edgevectorlength[edgeI]/edgelength[edgeI];
		directions.append(vector_dir*sign/std::abs(unitary_vector&vector_dir));
		heights.append(edgecentres[edgeI].z());
	};

	const labelUList edgeFaces = faPatchDict.edgeFaces(); 
	
	forAll(edgeFaces, edgeI)
	{
		vector nEdge = n[edgeFaces[edgeI]];
		slope.append(180.0/M_PI*Foam::atan(Foam::sqrt(sqr(nEdge.x())+sqr(nEdge.y()))/std::abs(nEdge.z())));
	};
	
	subHidrogram.add("wide", wide);
	subHidrogram.add("slope", slope);
	subHidrogram.add("directions", directions);
	subHidrogram.add("heights", heights);
	subHidrogram.add("alpha", sumwide/gHidrograms.subDict(boundaryName).lookupOrDefault<scalar>("wide",SMALL));
	subHidrogram.add("times", gHidrograms.subDict(boundaryName).lookup("times"));
	subHidrogram.add("flows", gHidrograms.subDict(boundaryName).lookup("flows"));
	subHidrogram.add("Cvs", gHidrograms.subDict(boundaryName).lookup("Cvs"));
	
	Hidrograms.add(name,subHidrogram);	
};

Info << "Hidrograms =" << Hidrograms << endl;
Info << "Creating initial flow condition " << endl;

#include "createFlowCondition.H"

h = max(h,hmin);
rho = max(rho,rho_w);
rhoh = rho*h;
Cw = (rho_s-rho)/(rho_s-rho_w);
Cw = min(dimensionedScalar("one", dimless, 1),max(dimensionedScalar("zero", dimless, 0),Cw));
Cwh = Cw*h;
Cv = 1-Cw;

h.write();
Cv.write();
Us.write();
Q.write();
c.write();
n.write();
A.write();

if (terrainModification)
{
	deltac0.write();
	he.write();
}
else
{
	deltah0.write();
}

if(Pstream::parRun())
{
    edgeVectorField ec
    (
        IOobject
        (
            "ec",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh.edgeCentres()
    );

    labelIOList edgeOwner
    (
        IOobject
        (
            "edgeOwner",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh.edgeOwner()
    );
	
    ec.write();
    edgeOwner.write();
}
