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
#include "Uchida.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace depositionModels
{
    defineTypeNameAndDebug(Uchida, 0);
    addToRunTimeSelectionTable(depositionModel, Uchida, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::depositionModels::Uchida::Uchida
(
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
    depositionModel(type(), depositionProperties, Us, h, he, rho, pb, tau, Cv),
    rhob_("rho_b", depositionProperties),
    rhow_("rho_w", depositionProperties),
    rhos_("rho_s", depositionProperties),
    phi_("phi", coeffDict_),
    deltad_("deltad", coeffDict_),
    d_("d", coeffDict_),
    hmin_("hmin", coeffDict_),
    g_("g", dimAcceleration, 9.81)
{
    Info<< "    " << rhow_ << nl
        << "    " << rhos_ << nl
        << "    " << phi_ << nl
        << "    " << deltad_ << nl
        << "    " << d_ << nl
        << "    " << hmin_ << nl << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::areaScalarField&
Foam::depositionModels::Uchida::Sd() const
{
    areaVectorField n = Us_.mesh().faceAreaNormals();

    vector ii (1,0,0), jj (0,1,0), kk (0,0,1) ;
    dimensionedVector i (ii), j (jj), k (kk);

    areaScalarField theta = Foam::atan(Foam::sqrt(sqr(n&i)+sqr(n&j))/mag(n&k));
    areaScalarField Ce = rhow_*Foam::tan(theta)/((rhos_-rhow_)*(Foam::tan(phi_)-Foam::tan(theta))) ;
    Ce = min(max(Ce, dimensionedScalar("zero", dimless,0)), dimensionedScalar("one", dimless,1.0)) ; 

    Sd_ = deltad_*(Cv_-Ce)/(Cvb_)*mag(Us_)*h_/d_*dZ_;
    Sd_ = max(min(Sd_, h_*Cv_/Us_.db().time().deltaT()/Cvb_),dimensionedScalar("zero", dimVelocity,0));

    forAll(Sd_, i)
    {
        if (mag(Us_.oldTime()[i]) < VSMALL)
        {
            Sd_[i] = Zero;
        }
        if (mag(h_[i]) <= hmin_.value())
        {
            Sd_[i] = Zero;
        }
    }

    return Sd_;
}

bool Foam::depositionModels::Uchida::read
(
    const dictionary& depositionProperties
)
{
    readDict(type(), depositionProperties);

    coeffDict_.readEntry("Cvb", Cvb_);
    coeffDict_.readEntry("phi", phi_);
    coeffDict_.readEntry("deltad", deltad_);
    coeffDict_.readEntry("d", d_);
    coeffDict_.readEntry("hmin", hmin_);

    return true;
}


// ************************************************************************* //
