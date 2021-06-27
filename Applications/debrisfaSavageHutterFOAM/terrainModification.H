deltavol = A*deltah;
areaScalarField nz = n & vector(0,0,-1);
deltah /= nz;

deltah.correctBoundaryConditions();

scalar maxdeltah = max(deltah.ref()).value();
scalar mindeltah = min(deltah.ref()).value();

Info << "min(deltah) = " << mindeltah << " and max(deltah) = "<< maxdeltah << endl;

#include "smoother.H"

if (deltahMin.value() < maxdeltah || -deltahMin.value() > mindeltah)
{	
	int nIter = 0;

	while (nIter < nTerrainIter && (deltahMin.value() < maxdeltah || -deltahMin.value() > mindeltah) )
	{
		nIter ++;

		Info<< "Performing terrain modification... iteration number = " << nIter << endl;
		
		#include "globalFields.H"

		scalarField meshptsdeltah(aMesh.nPoints(), zero());
		scalarField meshptsh(aMesh.nPoints(), zero());

		forAll(meshptsdeltah, pI)
		{
			labelList facesList = N_A_List[pointAddressing[pI]];

			forAll(facesList, fI)
			{
				scalar alpha_ij = gA[facesList[fI]]/N_A[pointAddressing[pI]];
				meshptsh[pI] += gh[facesList[fI]]*alpha_ij;
				meshptsdeltah[pI] += gdeltah[facesList[fI]]*alpha_ij;
			}
		} 

		forAll(meshPts, pI)
		{
			if (meshptsh[pI] > hMinterrain.value())
			{
				PointListPolyMesh[meshPts[pI]].z() += meshptsdeltah[pI]*alpha_relax;
			}		
		}

		tmp<scalarField> tvol = mesh.movePoints(PointListPolyMesh);
		scalarField& vol = tvol.ref();

		const fvPatchList& patches = mesh.boundary();

		forAll(patches,patchi)
		{
			if (patches[patchi].name() == patchName) 
			{
				Volsed.primitiveFieldRef() = patches[patchi].patchSlice(vol);
			}
		}

		if (aMesh.movePoints())
		{
			Info<< "Moving faMesh..." << endl;
		}

		A.ref() = aMesh.S();
		nz = aMesh.faceAreaNormals() & vector(0,0,-1);

		deltavol += Volsed;
		deltah = deltavol/(A+dimensionedScalar("A_0", dimArea, SMALL))/nz;

		deltah.correctBoundaryConditions();

		maxdeltah = max(deltah.ref()).value();
		mindeltah = min(deltah.ref()).value();

		Info << "New points after point motion min(deltah) = " << mindeltah << " and new max(deltah) = "<< maxdeltah << endl;

		#include "smoother.H"				
	}

	c = aMesh.areaCentres();
	n = aMesh.faceAreaNormals();
	gn = g & n;
	gs = g - gn*n;
}

deltah = deltah*nz;
deltah.correctBoundaryConditions();

