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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::entrainmentModel> Foam::entrainmentModel::New
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
{
    const word modelName(entrainmentProperties.get<word>("entrainmentModel"));

    Info<< nl << "Selecting entrainment model "
        << modelName << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelName);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown entrainmentModel " << modelName << nl << nl
            << "Valid types :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<entrainmentModel>
    (
        cstrIter()(entrainmentProperties, Us, h, he, rho, pb, tau, Cv)
    );
}


// ************************************************************************* //
