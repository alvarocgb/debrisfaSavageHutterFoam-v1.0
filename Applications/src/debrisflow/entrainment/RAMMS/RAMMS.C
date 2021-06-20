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

#include "fvCFD.H"
#include "faCFD.H"
#include "RAMMS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace entrainmentModels
{
    defineTypeNameAndDebug(RAMMS, 0);
    addToRunTimeSelectionTable(entrainmentModel,RAMMS, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::entrainmentModels::RAMMS::RAMMS
(
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
    entrainmentModel(type(), entrainmentProperties, Us, h, he, rho, pb, tau, Cv),
    dzdt_("dzdt", dimVelocity, coeffDict_),
    dzdtau_("dzdtau", dimAntimPa, coeffDict_),
    tauc_("tauc", dimPressure, coeffDict_),
    gs_(Us_.db().lookupObject<areaVectorField>("gs")),
    gn_(Us_.db().lookupObject<areaScalarField>("gn")),
    emm_
    (
        IOobject
        (
            "emm",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimLength)
    ),
    ett_
    (
        IOobject
        (
            "ett",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimLength)
    ),
    he0_
    (
        IOobject
        (
            "he0",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        he_
    )
{
    Info << "    " << dzdt_ <<  endl;
    Info << "    " << dzdtau_ << endl;
    Info << "    " << tauc_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::areaScalarField&
Foam::entrainmentModels::RAMMS::Sm() const
{
    ett_ = he0_-he_;

    Sm_ = dzdt_*eZ_;
    Sm_ = max(min(Sm_, he_/Us_.db().time().deltaT()),dimensionedScalar("zero", dimVelocity,0));

    forAll(Sm_, i)
    {    	
        if (mag(tau_[i]) <= tauc_.value())
        {
            Sm_[i] = Zero;
        }
	else
	{
	    emm_[i] = (mag(tau_[i])-tauc_.value())*dzdtau_.value(); 

	    if (ett_[i] > emm_[i])
	    {
		Sm_[i] = 0;
	    }
	}
    }
	
    return Sm_;
}


bool Foam::entrainmentModels::RAMMS::read
(
    const dictionary& entrainmentProperties
)
{
    readDict(type(), entrainmentProperties);

    coeffDict_.readEntry("dzdt", dzdt_);
    coeffDict_.readEntry("dzdtau", dzdtau_);
    coeffDict_.readEntry("tauc", tauc_);

    return true;
}


// ************************************************************************* //
