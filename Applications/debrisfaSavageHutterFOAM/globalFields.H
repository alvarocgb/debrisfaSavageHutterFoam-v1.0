if(Pstream::parRun())
{
	gh = Field<scalar>(gA.size(), zero());
	gdeltah = Field<scalar>(gA.size(), zero());

	forAll(faceAddressingown, fI)
	{
		gh[faceAddressingown[fI]] = h[fI];
		gdeltah[faceAddressingown[fI]] = deltah[fI];		
	}

	forAll (processorFAList, pI)
	{
		processorFaPatch pfaPatch
		(
			processorFAList[pI].get<word>("name"),
			processorFAList[pI],
			processorFAList[pI].get<label>("index"),
			faBoundaries
		);
	
		scalarList deltahown(deltah.ref().field());
		pfaPatch.send(Pstream::commsTypes::blocking,deltahown);
	}
	
	forAll (processorFAList, pI)
	{
		processorFaPatch pfaPatch
		(
			processorFAList[pI].get<word>("name"),
			processorFAList[pI],
			processorFAList[pI].get<label>("index"),
			faBoundaries
		);
	
		scalarList deltahngb = scalarList(processorFAList[pI].get<label>("nfacesngb"), zero());
		pfaPatch.receive(Pstream::commsTypes::blocking, deltahngb);

		labelList faceAddressingngb = processorFAList[pI].get<labelList>("faceAddressingngb");

		forAll(faceAddressingngb, fI)
		{
			gdeltah[faceAddressingngb[fI]] = deltahngb[fI];
		}
	}

	forAll (processorFAList, pI)
	{
		processorFaPatch pfaPatch
		(
			processorFAList[pI].get<word>("name"),
			processorFAList[pI],
			processorFAList[pI].get<label>("index"),
			faBoundaries
		);
	
		scalarList hown(h.ref().field());
		pfaPatch.send(Pstream::commsTypes::blocking,hown);
	}

	forAll (processorFAList, pI)
	{
		processorFaPatch pfaPatch
		(
			processorFAList[pI].get<word>("name"),
			processorFAList[pI],
			processorFAList[pI].get<label>("index"),
			faBoundaries
		);
	
		scalarList hngb = scalarList(processorFAList[pI].get<label>("nfacesngb"), zero());
		pfaPatch.receive(Pstream::commsTypes::blocking,hngb);

		labelList faceAddressingngb = processorFAList[pI].get<labelList>("faceAddressingngb");

		forAll(faceAddressingngb, fI)
		{
			gh[faceAddressingngb[fI]] = hngb[fI];
		}
	}		
}
else
{
	gh = h.ref().field();
	gdeltah = deltah.ref().field();
}
