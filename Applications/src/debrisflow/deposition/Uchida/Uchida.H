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

Class
    Foam::depositionModels::Uchida

Description
     

SourceFiles
    Uchida.C

Author
    

\*---------------------------------------------------------------------------*/

#ifndef Uchida_H
#define Uchida_H

#include "depositionModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"
#include "faMatrices.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace depositionModels
{

/*---------------------------------------------------------------------------*\
                           Class Uchida Declaration
\*---------------------------------------------------------------------------*/

class Uchida
:
    public depositionModel
{
    // Private data

        //- 
        dimensionedScalar hmin_;

        //- 
        dimensionedScalar phi_;

        //- 
        dimensionedScalar deltad_;

        //- 
        dimensionedScalar d_;

	//-
	dimensionedScalar rhow_;

	//-
	dimensionedScalar rhos_;

        //- 
        dimensionedScalar rhob_;

	//-
        dimensionedScalar g_;


public:

    //- Runtime type information
    TypeName("Uchida");


    // Constructors

        //- Construct from components
        Uchida
        (
            const dictionary& frictionProperties,
            const areaVectorField& Us,
            const areaScalarField& h,
            const areaScalarField& he,
            const areaScalarField& rho,
            const areaScalarField& pb,
            const areaVectorField& tau,
    	    const areaScalarField& Cv
        );


    //- Destructor
    virtual ~Uchida() = default;


    // Member Functions

        //- Return Sink by deposition
        virtual const areaScalarField& Sd() const;

        //- Read depositionProperties dictionary
        virtual bool read(const dictionary& depositionProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frictionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
