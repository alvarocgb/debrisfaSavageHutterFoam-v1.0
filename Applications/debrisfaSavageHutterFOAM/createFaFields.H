Info << "Reading field h" << endl;

areaScalarField h
(
	IOobject
	(
		"h",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	aMesh
);

Info << "Reading field Us" << endl;

areaVectorField Us
(
	IOobject
	(
		"Us",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	aMesh
);

areaScalarField Cv
(
	IOobject
	(
		"Cv",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	aMesh,
	dimless
);

areaScalarField Cw
(
	IOobject
	(
		"Cw",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	),
	aMesh,
	dimless
);

areaScalarField rho
(
	IOobject
	(
		"rho",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	),
	aMesh,
	dimDensity
);

areaScalarField rhoh
(
	IOobject
	(
		"rhoh",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	),
	aMesh
);

areaScalarField Cwh
(
	IOobject
	(
		"Cwh",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	),
	aMesh
);

edgeScalarField phi2s
(
	IOobject
	(
		"phi2s",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	fac::interpolate(rho*h*Us) & aMesh.Le()
);

edgeScalarField Q
(
	IOobject
	(
		"Q",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	fac::interpolate(h*Us) & aMesh.Le()
);

edgeScalarField phis
(
	IOobject
	(
		"phis",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	fac::interpolate(Us) & aMesh.Le()
);

areaScalarField A
(
	IOobject
	(
		"A",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	),
	aMesh
);

areaVectorField c
(
	IOobject
	(
		"c",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	aMesh.areaCentres()
);

areaVectorField n
(
	IOobject
	(
		"n",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	aMesh.faceAreaNormals()
);

areaScalarField gn
(
	IOobject
	(
		"gn",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	g & n
);

areaVectorField gs
(
	IOobject
	(
		"gs",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	g - gn*n
);

areaScalarField pb
(
	IOobject
	(
		"pb",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	aMesh,
	dimensionedScalar(dimPressure)
);

areaVectorField tau
(
	IOobject
	(
		"tau",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	aMesh,
	dimensionedVector(dimPressure)
);

Info << "Reading field he" << endl;

areaScalarField he
(
	IOobject
	(
		"he",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	),
	aMesh,
	dimensionedScalar(dimLength)
);

areaScalarField he0
(
	IOobject
	(
		"he0",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	he
);

areaScalarField deltah
(
	IOobject
	(
		"deltah",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	aMesh,
	dimensionedScalar(dimLength)
);

areaScalarField deltahh
(
	IOobject
	(
		"deltahh",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	aMesh,
	dimensionedScalar(dimLength)
);

areaScalarField deltah0
(
	IOobject
	(
		"deltah0",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE 
	),
	aMesh,
	dimensionedScalar(dimLength)
);

areaScalarField c0
(
	IOobject
	(
		"c0",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	aMesh,
	dimensionedScalar(dimLength)
);

c0 = c & vector(0,0,1);

areaScalarField deltac0
(
	IOobject
	(
		"deltac0",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	aMesh,
	dimensionedScalar(dimLength)
);
