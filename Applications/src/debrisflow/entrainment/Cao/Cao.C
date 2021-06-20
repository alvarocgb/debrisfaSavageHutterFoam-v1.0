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
#include "Cao.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace entrainmentModels
{
    defineTypeNameAndDebug(Cao, 0);
    addToRunTimeSelectionTable(entrainmentModel, Cao, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::entrainmentModels::Cao::Cao
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
    h0_("h0", entrainmentProperties),
    d_("d", coeffDict_),
    beta_("beta", coeffDict_),
    tauc_("tauc", dimPressure, coeffDict_),
    g_("g", dimAcceleration, 9.81),
    gs_(Us_.db().lookupObject<areaVectorField>("gs")),
    gn_(Us_.db().lookupObject<areaScalarField>("gn"))
{
    Info << "    " << d_ <<  endl;
    Info << "    " << h0_ <<  endl;
    Info << "    " << beta_ <<  endl;
    Info << "    " << tauc_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::areaScalarField&
Foam::entrainmentModels::Cao::Sm() const
{
    dimensionedScalar s = (rhos_-rhow_)/rhow_;

    areaScalarField theta = mag(tau_)/(rho_*g_*s*d_);
    areaScalarField thetac = mag(tauc_)/(rho_*g_*s*d_);

    Sm_ = beta_*(theta-thetac)*mag(Us_)/(h_+h0_)*pow(d_, -0.2)*eZ_;
    Sm_ = max(min(Sm_, he_/Us_.db().time().deltaT()),dimensionedScalar("zero", dimVelocity,0));

    return Sm_;
}


bool Foam::entrainmentModels::Cao::read
(
    const dictionary& entrainmentProperties
)
{
    readDict(type(), entrainmentProperties);

    coeffDict_.readEntry("d", d_);
    coeffDict_.readEntry("beta", beta_);
    coeffDict_.readEntry("tauc", tauc_);

    return true;
}


// ************************************************************************* //
