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
    Foam::frictionModels::VoellmyManning

Description
    

SourceFiles
    VoellmyManning.C

Author

\*---------------------------------------------------------------------------*/

#ifndef VoellmyManning_H
#define VoellmyManning_H

#include "frictionModel.H"
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
namespace frictionModels
{

/*---------------------------------------------------------------------------*\
                           Class VoellmyManning Declaration
\*---------------------------------------------------------------------------*/

const dimensionSet dimManning(0,-0.333333,1, 0, 0, 0, 0);
const dimensionSet dimAntiL(0,-1,0, 0, 0, 0, 0);

class VoellmyManning
:
    public frictionModel
{
    // Private data

        //- Coulomb friction parameter mu
        areaScalarField mu_;

        //- Voellmy friction parameter xi
        areaScalarField xi_;

        //- Equivalent Manning's n-value
        areaScalarField n_;

        //- Gravitational acceleration
    	dimensionedScalar g_;

        //- 
        areaScalarField Cvlim_;

public:

    //- Runtime type information
    TypeName("VoellmyManning");


    // Constructors

        //- Construct from components
        VoellmyManning
        (
            const dictionary& frictionProperties,
            const areaVectorField& Us,
            const areaScalarField& h,
            const areaScalarField& rho,
            const areaScalarField& p,
    	    const areaScalarField& Cv
        );


    //- Destructor
    virtual ~VoellmyManning() = default;


    // Member Functions

        //- Return explicit coefficients for the linear system of equations
        //- considering friction
        virtual const areaVectorField& tauSc() const;

        //- Return implicit coefficients for the linear system of equations
        //- considering friction
        virtual const areaScalarField& tauSp() const;

        //- Read the friction properties dictionary
        virtual bool read(const dictionary& frictionProperties);

        //- Return the flow height used as boundary condition for the boundary patch
        virtual const List<scalar>& ht(const float flux, const List<scalar> wide, const List<scalar> slope, const List<scalar> height, const float rhot, const word bdName) const;

        //- Return the uniform flow height used as boundary condition for the boundary patch
        virtual const scalar& htu(const float flux, const List<scalar> wide, const List<scalar> slope, const float rhot, const word bdName) const;

	//- Return the flow velocity used as boundary condition for the boundary patch
        virtual const scalar& ut(const float flux, const List<scalar> wide, const List<scalar> slope, const float rhot, const word bdName) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frictionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
