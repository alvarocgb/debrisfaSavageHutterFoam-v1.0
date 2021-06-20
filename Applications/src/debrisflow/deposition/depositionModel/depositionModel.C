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

#include "depositionModel.H"
#include "fvCFD.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(depositionModel, 0);
    defineRunTimeSelectionTable(depositionModel, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::depositionModel::readDict
(
    const word& type,
    const dictionary& dict
)
{
    depositionProperties_ = dict;
    coeffDict_ = depositionProperties_.optionalSubDict(type + "Coeffs");
}


void Foam::depositionModel::readdepositionZones()
{
    word deposition = coeffDict_.lookupOrDefault<word>("deposition", "uniform");
   
    Info << "type " << deposition << endl;

    if (deposition == "uniform")
    {
        dZ_ = dimensionedScalar("uniform", dimless, 1.0);
    }
    else if (deposition == "nonuniform")
    {
        Info<< "Reading deposition zones" << endl;

	areaVectorField c = Us_.mesh().areaCentres();

        const PtrList<entry> regions
        (
            coeffDict_.lookup("depositionZones")
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
                	dZ_[i] = 1;
                }
            }        
        }	
    }
    dZ_.write();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::depositionModel::depositionModel
(
    const word& type,
    const dictionary& depositionProperties,
    const areaVectorField& Us,
    const areaScalarField& h,
    const areaScalarField& he,
    const areaScalarField& rho,
    const areaScalarField& pb,
    const areaVectorField& tau,
    const areaScalarField& Cv
)
:
    depositionProperties_(depositionProperties),
    coeffDict_
    (
        depositionProperties_.optionalSubDict(type + "Coeffs")
    ),
    rho_(rho),
    Us_(Us),
    h_(h),
    he_(he),
    pb_(pb),
    tau_(tau),
    Cv_(Cv),
    Cvb_(dimensionedScalar("Cvb", dimless,0)),
    Sd_
    (
        IOobject
        (
            "Sd",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimVelocity)
    ),
    dZ_
    (
        IOobject
        (
            "dZ",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimless)
    )
{
    readdepositionZones();
    dimensionedScalar rho_w("rho_w", dimDensity, depositionProperties);
    dimensionedScalar rho_s("rho_s", dimDensity, depositionProperties);
    dimensionedScalar rho_b("rho_b", dimDensity, depositionProperties);
    Cvb_ = dimensionedScalar("Cvb", dimless, (rho_b.value()-rho_w.value())/(rho_s.value()-rho_w.value()));

    Info << "    " << Cvb_ << endl;
}


// ************************************************************************* //
