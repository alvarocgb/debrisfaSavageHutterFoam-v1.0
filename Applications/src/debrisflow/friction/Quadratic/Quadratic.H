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
    Foam::frictionModels::Quadratic

Description
    An incompressible Quadratic friction model like the one used in FLO-2D

SourceFiles
    Quadratic.C

Author
    Álvaro González Bilbao

\*---------------------------------------------------------------------------*/

#ifndef Quadratic_H
#define Quadratic_H

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
                           Class Quadratic Declaration
\*---------------------------------------------------------------------------*/

const dimensionSet dimManning(0,-0.333333,1, 0, 0, 0, 0);
const dimensionSet dimAntiL(0,-1,0, 0, 0, 0, 0);

class Quadratic
:
    public frictionModel
{
    // Private data


        //- Resistance parameter K
        areaScalarField K_;

        //- Equivalent Manning's n-value
        areaScalarField n_;

        //-  
        areaScalarField alpha1_;

        //- 
        areaScalarField alpha2_;

        //-  
        areaScalarField beta1_;

        //- 
        areaScalarField beta2_;

        //- Gravitational acceleration
    	dimensionedScalar g_;

        //- Dynamic viscosity of the fluid matrix
        mutable areaScalarField eta_;

        //- Yield stress of the fluid matrix (tau_y)
        mutable areaScalarField Ty_;

    // Protected Member Functions

        //- Set Ty field to zero
        void resetTy() const;

        //- Set eta field to zero
        void reseteta() const;

        //- Calculate Ty field
        void calcTy() const;

        //- Calculate eta field
        void calceta() const;

public:

    //- Runtime type information
    TypeName("Quadratic");


    // Constructors

        //- Construct from components
        Quadratic
        (
            const dictionary& frictionProperties,
            const areaVectorField& Us,
            const areaScalarField& h,
            const areaScalarField& rho,
            const areaScalarField& p,
            const areaScalarField& Cv
        );


    //- Destructor
    virtual ~Quadratic() = default;


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
