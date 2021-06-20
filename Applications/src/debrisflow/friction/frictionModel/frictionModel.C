/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | debrisfaSavageHutterFOAM
    \\  /    A nd           | Copyright (C)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Author


\*---------------------------------------------------------------------------*/

#include "frictionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(frictionModel, 0);
    defineRunTimeSelectionTable(frictionModel, dictionary);
}

// file-local
static const Foam::dimensionSet dimTauSc(1, -1, -2, 0, 0, 0, 0);
static const Foam::dimensionSet dimTauSp(1, -2, -1, 0, 0, 0, 0);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::frictionModel::readDict
(
    const word& type,
    const dictionary& dict
)
{
    frictionProperties_ = dict;
    coeffDict_ = frictionProperties_.optionalSubDict(type + "Coeffs");
}


void Foam::frictionModel::resetTauSc() const
{
    tauSc_ = dimensionedVector(dimTauSc);
}


void Foam::frictionModel::resetTauSp() const
{
    tauSp_ = dimensionedScalar(dimTauSp);
}


void Foam::frictionModel::readfield
(
    areaScalarField& fp, 
    word fpName
)
{
    const PtrList<entry> regions
    (
        coeffDict_.lookup("regions")
    );

    wordList areaNames;
    areaNames.setSize(regions.size());

    Info << nl << "reading " + fpName + " ..." << endl;

    forAll(regions, areaI)
    {
	areaVectorField c = Us_.mesh().areaCentres();

        const entry& regionInfo = regions[areaI];

        areaNames[areaI] = regionInfo.keyword();

        dictionary areaDict = regionInfo.dict();

	word type;
        dimensionedScalar fpi;	

        areaDict.lookup(fpName) >> fpi;
	areaDict.lookup("type") >> type;

	Info << fpName << "  value in area "<< areaNames[areaI] 
	     << " is " << fpi.value() <<endl;

	if (type == "default")
	{
            forAll(c.internalField(), i)
            {
        	fp[i] = fpi.value();
            }
	}
	else if (type == "polygon")
	{
            vector offset;
            List<point2D> points;
            List<vector> vertices;

            areaDict.lookup("offset") >> offset;
            areaDict.lookup("vertices") >> vertices;

            points.resize(vertices.size());

            forAll(vertices, vI)
            {
    	        points[vI] = point2D(vertices[vI].x()+offset.x(), vertices[vI].y()+offset.y());
            }

            HormannAgathos polygon(points, 0.001);

            forAll(c.internalField(), i)
            {
    	        if (polygon.evaluate(point2D(c[i].x(), c[i].y())) != HormannAgathos::POINT_OUTSIDE)
                {
        	        fp[i] = fpi.value();
                }
            }
	}
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::frictionModel::frictionModel
(
    const word& type,
    const dictionary& frictionProperties,
    const areaVectorField& Us,
    const areaScalarField& h,
    const areaScalarField& rho,
    const areaScalarField& p,
    const areaScalarField& Cv
)
:
    frictionProperties_(frictionProperties),
    coeffDict_
    (
        frictionProperties_.optionalSubDict(type + "Coeffs")
    ),
    u0_("u0", dimVelocity, frictionProperties_),
    h0_("h0", dimLength, frictionProperties_),
    Us_(Us),
    h_(h),
    rho_(rho),
    p_(p),
    Cv_(Cv),
    tauSp_
    (
        IOobject
        (
            "tauSp",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimTauSp)
    ),
    tauSc_
    (
        IOobject
        (
            "tauSc",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedVector(dimTauSc)
    )
{
    Info<< "    with " << nl
        << "    " << u0_ << nl
        << "    " << h0_ << endl;
}


// ************************************************************************* //
