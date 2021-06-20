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

#include "entrainmentModel.H"
#include "fvCFD.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(entrainmentModel, 0);
    defineRunTimeSelectionTable(entrainmentModel, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::entrainmentModel::readDict
(
    const word& type,
    const dictionary& dict
)
{
    entrainmentProperties_ = dict;
    coeffDict_ = entrainmentProperties_.optionalSubDict(type + "Coeffs");
}

void Foam::entrainmentModel::readentrainmentZones()
{
    word entrainment = coeffDict_.lookupOrDefault<word>("entrainment", "uniform");

    if (entrainment == "uniform")
    {
        eZ_ = dimensionedScalar("uniform", dimless, 1.0);
    }
    else if (entrainment == "nonuniform")
    {
        Info<< "Reading entrainment zones" << endl;

	areaVectorField c = Us_.mesh().areaCentres();

        const PtrList<entry> regions
        (
            coeffDict_.lookup("entrainmentZones")
        );

        wordList areaNames;
        areaNames.setSize(regions.size());

        forAll(regions, areaI)
        {
            const entry& regionInfo = regions[areaI];

            Info<< "processing region " << regions[areaI].keyword() << endl;

            areaNames[areaI] = regionInfo.keyword();

            dictionary areaDict = regionInfo.dict();

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
                	eZ_[i] = 1;
                }
            }        
        }	
    }
    eZ_.write();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::entrainmentModel::entrainmentModel
(
    const word& type,
    const dictionary& entrainmentProperties,
    const areaVectorField& Us,
    const areaScalarField& h,
    const areaScalarField& he,
    const areaScalarField& rho,
    const areaScalarField& pb,
    const areaVectorField& tau,
    const areaScalarField& Cv
)
:
    entrainmentProperties_(entrainmentProperties),
    coeffDict_
    (
        entrainmentProperties_.optionalSubDict(type + "Coeffs")
    ),
    rho_(rho),
    Cv_(Cv),
    Us_(Us),
    h_(h),
    he_(he),
    pb_(pb),
    tau_(tau),
    Sm_
    (
        IOobject
        (
            "Sm",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimVelocity)
    ),
    eZ_
    (
        IOobject
        (
            "eZ",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimless)
    )
{
    readentrainmentZones();
}


// ************************************************************************* //
