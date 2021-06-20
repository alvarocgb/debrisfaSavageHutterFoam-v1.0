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
#include "Egashira.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace entrainmentModels
{
    defineTypeNameAndDebug(Egashira, 0);
    addToRunTimeSelectionTable(entrainmentModel, Egashira, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::entrainmentModels::Egashira::Egashira
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
    rhow_("rho_w", entrainmentProperties),
    rhos_("rho_s", entrainmentProperties),
    phi_("phi", coeffDict_),
    tauc_("tauc", dimPressure, coeffDict_),
    gs_(Us_.db().lookupObject<areaVectorField>("gs")),
    gn_(Us_.db().lookupObject<areaScalarField>("gn"))
{
    Info << "    " << phi_ <<  endl;
    Info << "    " << tauc_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::areaScalarField&
Foam::entrainmentModels::Egashira::Sm() const
{
    areaVectorField n = Us_.mesh().faceAreaNormals();

    areaScalarField thetae = Foam::atan(Cv_*(rhos_-rhow_)/(Cv_*(rhos_-rhow_)+rhow_)*Foam::tan(M_PI*phi_/180));

    vector ii (1,0,0), jj (0,1,0), kk (0,0,1) ;

    dimensionedVector i (ii), j (jj), k (kk);

    areaScalarField theta = Foam::atan(Foam::sqrt(sqr(n&i)+sqr(n&j))/mag(n&k));

    Sm_ = mag(Us_)*Foam::tan(theta-thetae)*eZ_;
    Sm_ = max(min(Sm_, he_/Us_.db().time().deltaT()),dimensionedScalar("zero", dimVelocity,0));

    forAll(Sm_, i)
    {
        if (mag(tau_[i]) < tauc_.value())
        {
            Sm_[i] = Zero;
        }
    }


    return Sm_;
}


bool Foam::entrainmentModels::Egashira::read
(
    const dictionary& entrainmentProperties
)
{
    readDict(type(), entrainmentProperties);

    coeffDict_.readEntry("phi", phi_);
    coeffDict_.readEntry("tauc", tauc_);

    return true;
}


// ************************************************************************* //
