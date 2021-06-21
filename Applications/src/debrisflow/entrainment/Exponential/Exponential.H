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
    Foam::entrainmentModels::Exponential

Description
    

SourceFiles
    Exponential.C

Author
    

\*---------------------------------------------------------------------------*/

#ifndef Exponential_H
#define Exponential_H

#include "entrainmentModel.H"
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
namespace entrainmentModels
{

/*---------------------------------------------------------------------------*\
                           Class Exponential Declaration
\*---------------------------------------------------------------------------*/
const dimensionSet dimAntiL(0,-1,0, 0, 0, 0, 0);

class Exponential
:
    public entrainmentModel
{
    // Private data

        //-
        dimensionedScalar E_;

        //- 
        dimensionedScalar tauc_;	

        //- Reference to the gravitation fields
        const areaVectorField& gs_;
        const areaScalarField& gn_;

public:

    //- Runtime type information
    TypeName("Exponential");


    // Constructors

        //- Construct from components
        Exponential
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
    virtual ~Exponential() = default;


    // Member Functions

        //- Return the Source by entrainment
        virtual const areaScalarField& Sm() const;

        //- Read entrainmentProperties dictionary
        virtual bool read(const dictionary& entrainmentProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frictionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
