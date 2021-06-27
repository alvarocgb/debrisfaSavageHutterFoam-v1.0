if (deltahMin.value()*alpha_smooth < maxdeltah || -deltahMin.value()*alpha_smooth > mindeltah)
{
	Info << "Smoothing deltah field..." << endl;

	//areaScalarField deltah_0 = deltah;

	List<List<label>> M_faces (aMesh.nFaces());
	List<List<label>> M_edges (aMesh.nFaces());

	scalarList deltahedge = List<scalar>();
	//scalarList heedge = List<scalar>();
	scalarList Aedge = List<scalar>();
	
	if(Pstream::parRun())
	{
		deltahedge = List<scalar>(aMesh.nEdges()-aMesh.nInternalEdges(), zero());
		areaScalarField::Boundary boundarydeltah(deltah.boundaryFieldRef());

		//heedge = List<scalar>(aMesh.nEdges()-aMesh.nInternalEdges(), zero());
		//areaScalarField::Boundary boundaryhe(he.boundaryFieldRef());

		Aedge = List<scalar>(aMesh.nEdges()-aMesh.nInternalEdges(), zero());
		areaScalarField::Boundary boundaryA(A.boundaryFieldRef());

		forAll(faboundaryNames, faboundaryI)
		{		
			if (faboundaryTypes[faboundaryI] == "processor")
			{
				labelList bdList(faBoundaries[faboundaryI]);
				List<scalar> bdfielddeltah(boundarydeltah[faboundaryI]);
				//List<scalar> bdfieldhe(boundaryhe[faboundaryI]);
				List<scalar> bdfieldA(boundaryA[faboundaryI]);				

				forAll(bdList, edgeI)
				{
					deltahedge[bdList[edgeI]-aMesh.nInternalEdges()] = bdfielddeltah[edgeI];
					//heedge[bdList[edgeI]-aMesh.nInternalEdges()] = bdfieldhe[edgeI];
					Aedge[bdList[edgeI]-aMesh.nInternalEdges()] = bdfieldA[edgeI];
				}
			}
		}
	}

	forAll (faceNeighbours, faceI)
	{
		List<label> faceList = faceNeighbours[faceI];
		List<label> edgeList = edgeNeighbours[faceI];

		scalar min = deltah[faceI];
		scalar max = deltah[faceI];

		averagedeltah[faceI] = deltah[faceI]/(NfaceNeighbours[faceI]+NedgeNeighbours[faceI]+1);
		
		forAll(faceList, fI)
		{
			averagedeltah[faceI] += deltah[faceList[fI]]/(NfaceNeighbours[faceI]+NedgeNeighbours[faceI]+1);

			if (deltah[faceList[fI]]< min) { min =  deltah[faceList[fI]]; }
			if (deltah[faceList[fI]]> max) { max =  deltah[faceList[fI]]; }		
		}

		forAll(edgeList, eI)
		{
			averagedeltah[faceI] += deltahedge[edgeList[eI]-aMesh.nInternalEdges()]/(NfaceNeighbours[faceI]+NedgeNeighbours[faceI]+1);

			if (deltahedge[edgeList[eI]-aMesh.nInternalEdges()]< min) { min = deltahedge[edgeList[eI]-aMesh.nInternalEdges()]; }
			if (deltahedge[edgeList[eI]-aMesh.nInternalEdges()]> max) { max = deltahedge[edgeList[eI]-aMesh.nInternalEdges()]; }		
		}

		sigmadeltah[faceI] = max-min;

		if (std::abs(deltah[faceI]-averagedeltah[faceI]) <= alpha_gamma*sigmadeltah[faceI])
		{
			Gamma[faceI] = 0;
		}
/*		else if ((deltah[faceI] > averagedeltah[faceI]) && (he[faceI] <= 0))
		{
			Gamma[faceI] = 0;
		}*/
		else
		{
			Gamma[faceI] = 1;

			if (deltah[faceI] > averagedeltah[faceI])
			{
				forAll(faceList, fI)
				{
					if (deltah[faceI] > deltah[faceList[fI]])
					{
						M_faces[faceI].append(faceList[fI]);
					}
				}

				forAll(edgeList, eI)
				{
					if (deltah[faceI] > deltahedge[edgeList[eI]-aMesh.nInternalEdges()])
					{
						M_edges[faceI].append(edgeList[eI]);
					}
				}
			}
			else if (deltah[faceI] < averagedeltah[faceI])
			{
				forAll(faceList, fI)
				{
					if (deltah[faceI] < deltah[faceList[fI]]) //&& (he[faceList[fI]] > 0))
					{
						M_faces[faceI].append(faceList[fI]);
					}
				}

				forAll(edgeList, eI)
				{
					if (deltah[faceI] < deltahedge[edgeList[eI]-aMesh.nInternalEdges()]) //&& (heedge[edgeList[eI]-aMesh.nInternalEdges()] > 0))
					{
						M_edges[faceI].append(edgeList[eI]);
					}
				}
			}
		}
	}

	deltaVolMax = (deltah-averagedeltah)*beta_gamma*Gamma*A*nz;

	forAll (faceNeighbours, faceI)
	{
		if (Gamma[faceI] == 1)
		{
			List<label> MList_faces = M_faces[faceI];
			List<label> faceList = faceNeighbours[faceI];
			List<label> edgeList = edgeNeighbours[faceI];
		
			scalar sum_AN = 0;	
		
			forAll(faceList, fI)
			{
				sum_AN += A[faceList[fI]];
			}

			forAll(edgeList, eI)
			{
				sum_AN += Aedge[edgeList[eI]-aMesh.nInternalEdges()];
			}

			forAll(MList_faces, fI)
			{
				scalar alphaij = A[MList_faces[fI]]/sum_AN;
				deltavol[MList_faces[fI]] += deltaVolMax[faceI]*alphaij;
				deltavol[faceI] -= deltaVolMax[faceI]*alphaij;
			}	
		}
	}

	if(Pstream::parRun())
	{
		forAll(processorFAList, pI)
		{
			processorFaPatch pfaPatch
			(
				processorFAList[pI].get<word>("name"),
				processorFAList[pI],
				processorFAList[pI].get<label>("index"),
				faBoundaries
			);

			labelUList edgeFaces = pfaPatch.edgeFaces();
			labelList edgeList = processorFAList[pI].get<List<label>>("edgeLabels");

			scalarList deltavoledge = List<scalar>(edgeList.size(), zero());	
		
			forAll(edgeFaces, eI)
			{
				List<label> MList_edges = M_edges[edgeFaces[eI]];

				if (MList_edges.found(edgeList[eI]))
				{
					List<label> faceNList = faceNeighbours[edgeFaces[eI]];
					List<label> edgeNList = edgeNeighbours[edgeFaces[eI]];
				
					scalar sum_AN = 0;	
				
					forAll(faceNList, faceI)
					{
						sum_AN += A[faceNList[faceI]];
					}

					forAll(edgeNList, edgeI)
					{
						sum_AN += Aedge[edgeNList[edgeI]-aMesh.nInternalEdges()];
					}

					scalar alphaij = Aedge[edgeList[eI]-aMesh.nInternalEdges()]/sum_AN;
					deltavoledge[eI] += deltaVolMax[edgeFaces[eI]]*alphaij;
					deltavol[edgeFaces[eI]] -= deltaVolMax[edgeFaces[eI]]*alphaij;				
				}		
			}

			pfaPatch.send(Pstream::commsTypes::blocking,deltavoledge);
		}

		forAll(processorFAList, pI)
		{
			processorFaPatch pfaPatch
			(
				processorFAList[pI].get<word>("name"),
				processorFAList[pI],
				processorFAList[pI].get<label>("index"),
				faBoundaries
			);

			labelUList edgeFaces = pfaPatch.edgeFaces();

			scalarList deltavoledge = List<scalar>(edgeFaces.size(), zero());	
			pfaPatch.receive(Pstream::commsTypes::blocking, deltavoledge);

			forAll(edgeFaces, eI)
			{
				deltavol[edgeFaces[eI]] += deltavoledge[eI];
			}
		}
	}

	deltah = deltavol/(A+dimensionedScalar("A_0", dimArea, SMALL))/nz;

	//he = he + (deltah - deltah_0)*nz; //here there might be a problem, since there is no limit for the he value and then it can take negative values.

	maxdeltah = max(deltah.ref()).value();
	mindeltah = min(deltah.ref()).value();

	Info << "New values after smoothing are min(deltah) = " << mindeltah << " and max(deltah) = " << maxdeltah << endl;
}
