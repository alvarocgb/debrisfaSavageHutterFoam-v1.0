areaScalarField deltavol
(
	IOobject
	(
	    "deltavol",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	aMesh,
	dimensionedScalar(dimVolume)
);

areaScalarField Volsed
(
	IOobject
	(
	    "Volsed",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	aMesh,
	dimensionedScalar(dimVolume)
);

areaScalarField averagedeltah
(
	IOobject
	(
	    "averagedeltah",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	aMesh,
	dimensionedScalar(dimLength)
);

areaScalarField sigmadeltah
(
	IOobject
	(
	    "sigmadeltah",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	aMesh,
	dimensionedScalar(dimLength)
);

areaScalarField Gamma
(
	IOobject
	(
	    "Gamma",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	aMesh,
	dimensionedScalar(dimless)
);

areaScalarField deltaVolMax  
(
	IOobject
	(
	    "deltaVolMax",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	aMesh,
	dimensionedScalar(dimVolume)
);

pointField PointListPolyMesh = Field<vector>();
labelList meshPts = List<label>();
scalarField NfaceNeighbours = Field<scalar>();
List<List<label>> faceNeighbours = List<List<label>>();
scalarField NedgeNeighbours = Field<scalar>();
List<List<label>> edgeNeighbours = List<List<label>>();

scalarList N_A = List<scalar>();
List<List<label>> N_A_List = List<List<label>>();

scalarList gA = List<scalar>(); 
scalarField gh = Field<scalar>();
scalarField gdeltah = Field<scalar>();
labelList pointAddressing = List<label>();
labelList faceAddressingown = List<label>();
word patchName = word();
List<dictionary> processorFAList;

if (terrainModification)
{  
	Info << nl << "Creating Terrain Modification Fields" << endl;

	IOdictionary faMeshDefinition
	(
		IOobject
		(
			"faMeshDefinition",
			runTime.constant(),
			"faMesh",
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	wordList patchNameList = faMeshDefinition.lookup("polyMeshPatches");
	patchName = patchNameList[0];
	
	PointListPolyMesh = aMesh.mesh().points(); //raw points of the polyMesh object linked to the faMesh object
	meshPts = aMesh.patch().meshPoints(); //patch() returns constant reference to primitive patch. meshPoints return labelList of mesh points in patch

    	faceNeighbours  = List<List<label>>(aMesh.nFaces());
	NfaceNeighbours = Field<scalar>(aMesh.nFaces(), zero());
	edgeNeighbours  = List<List<label>>(aMesh.nFaces());
	NedgeNeighbours = Field<scalar>(aMesh.nFaces(), zero());

	IOdictionary N_A_dict //This is created by newreleaseAreaMapping when the terrainModification button is "on"
	(
	    IOobject
	    (
		"N_A_dict",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    )
	);

	N_A = N_A_dict.lookup("N_A"); //It is the number of faces to which every point belongs to.
	N_A_List = N_A_dict.lookupOrDefault<List<List<label>>>("N_A_List", List<List<label>>()); //Labels of the faces to which every point belongs to
	gA = N_A_dict.lookup("A"); //Areas of all the faces elements.	
	gh = Field<scalar>(gA.size(), zero());
	gdeltah = Field<scalar>(gA.size(), zero());

	if(Pstream::parRun())
	{
		labelIOList pointProcAddressing
		(
			IOobject
			(
				"pointProcAddressing",
				runTime.constant(),
				"faMesh",
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			)
		);
		pointAddressing = pointProcAddressing;
	}
	else
	{
		pointAddressing = labelList(aMesh.nPoints(), zero());
		forAll(pointAddressing, pI)
		{
			pointAddressing[pI] = pI;
		}
	}

	if(Pstream::parRun())
	{
		forAll(faboundaryNames, faboundaryI)
		{		
			if (faboundaryTypes[faboundaryI] == "processor")
			{
				labelList edgeList(faBoundaries[faboundaryI]);

				label own = 0;
				label ngb = 0;
				word name = faboundaryNames[faboundaryI];
				name.removeStart("procBoundary");
				read(name.substr(0,name.find("to")), own); 
				name.removeStart(name.substr(0,name.find("to"))+"to");			
				read(name, ngb);
				dictionary dictpatch(name);
				
				dictpatch.add("type", "processor");
				dictpatch.add("myProcNo", own);
				dictpatch.add("neighbProcNo", ngb);
				dictpatch.add("edgeLabels", edgeList);
				dictpatch.add("ngbPolyPatchIndex", faBoundaries[faboundaryI].ngbPolyPatchIndex());
				dictpatch.add("name", faboundaryNames[faboundaryI]);
				dictpatch.add("index", faboundaryI);

				processorFaPatch pfaPatch
				(
					faboundaryNames[faboundaryI],
					dictpatch,
					faboundaryI,
					faBoundaries
				);

				List<label> nfacesown = List<label>();
				nfacesown.append(aMesh.nFaces());
				pfaPatch.send(Pstream::commsTypes::blocking, nfacesown);

				List<label> nfacesngb = List<label>(1, zero());
				pfaPatch.receive(Pstream::commsTypes::blocking, nfacesngb);
				dictpatch.add("nfacesngb", nfacesngb[0]);

				labelIOList faceProcAddressing
				(
				    IOobject
				    (
					"faceProcAddressing",
					runTime.constant(),
					"faMesh",
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				    )
				);

				faceAddressingown = faceProcAddressing;
				pfaPatch.send(Pstream::commsTypes::blocking, faceAddressingown);

				labelList faceAddressingngb = labelList(nfacesngb[0], zero());
				pfaPatch.receive(Pstream::commsTypes::blocking, faceAddressingngb);
				dictpatch.add("faceAddressingngb", faceAddressingngb);

				processorFAList.append(dictpatch);
			}
		}
	}
	
	labelList edgeOwn = aMesh.edgeOwner();
	labelList edgeNgb = aMesh.edgeNeighbour();

	forAll(edgeNgb, edgeI)
	{
		label face_own = edgeOwn[edgeI];
		label face_ngb = edgeNgb[edgeI];
	
		NfaceNeighbours[face_own] ++;
		NfaceNeighbours[face_ngb] ++;

		faceNeighbours[face_own].append(face_ngb);
		faceNeighbours[face_ngb].append(face_own);
	}

	if(Pstream::parRun())
	{
		forAll(faboundaryNames, faboundaryI)
		{		
			if (faboundaryTypes[faboundaryI] == "processor")
			{
				labelList edgeList(faBoundaries[faboundaryI]);
				labelUList edgeFaces(faBoundaries[faboundaryI].edgeFaces());

				forAll(edgeList, edgeI)
				{
					NedgeNeighbours[edgeFaces[edgeI]] ++;
		    			edgeNeighbours[edgeFaces[edgeI]].append(edgeList[edgeI]);
				}		
			}
		}
	}
}
