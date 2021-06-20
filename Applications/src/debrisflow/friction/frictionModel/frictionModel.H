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

Namespace
    Foam::frictionModels

Description
    A namespace for various friction model implementations.

Class
    Foam::frictionModel

Description
    An abstract base class for friction models

SourceFiles
    frictionModel.C
    frictionModelNew.C

Author

\*---------------------------------------------------------------------------*/

#ifndef frictionModel_H
#define frictionModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "dimensionedScalar.H"
#include "tmp.H"
#include "autoPtr.H"
#include "faMatrices.H"
#include "areaFieldsFwd.H"
#include "FieldFields.H"
#include "HormannAgathos.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class frictionModel Declaration
\*---------------------------------------------------------------------------*/

class frictionModel
{
protected:

    // Protected data

        dictionary frictionProperties_;

        //- Model coefficients dictionary
        dictionary coeffDict_;

        //- Tolerance for velocity
        dimensionedScalar u0_;

        //- Tolerance for flow height
        dimensionedScalar h0_;

        //- Reference to the surface velocity field
        const areaVectorField& Us_;

        //- Reference to the flow height field
        const areaScalarField& h_;

        //- Reference to the flow density field
        const areaScalarField& rho_;

        //- Reference to the bottom pressure field
        const areaScalarField& p_;

        //- Reference to the volumetric sediment concentration
        const areaScalarField& Cv_;

        //- Implicit part of friction
        mutable areaScalarField tauSp_;

        //- Explicit part of friction
        mutable areaVectorField tauSc_;

        //- Calculated flow velocity used as boundary condition
        mutable scalar ut_;

        //- Calculated flow height used as boundary condition
        mutable List<scalar> ht_;


    // Protected Member Functions

        //- Read/update the frictionProperties and coeffDict dictionaries
        void readDict(const word& type, const dictionary& dict);

        //- Set tauSc field to zero
        void resetTauSc() const;

        //- Set tauSp field to zero
        void resetTauSp() const;


        //- Disallow copy construct
        frictionModel(const frictionModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const frictionModel&) = delete;


public:

    //- Runtime type information
    TypeName("frictionModel");


    // Declare run-time constructor selection table

#ifndef SWIG
        declareRunTimeSelectionTable
        (
            autoPtr,
            frictionModel,
            dictionary,
            (
                const dictionary& frictionProperties,
                const areaVectorField& Us,
                const areaScalarField& h,
                const areaScalarField& rho,
                const areaScalarField& p,
                const areaScalarField& Cv
            ),
            (frictionProperties, Us, h, rho, p, Cv)
        );
#endif


    // Selectors

        //- Return a reference to the selected friction model
        static autoPtr<frictionModel> New
        (
            const dictionary& frictionProperties,
            const areaVectorField& Us,
            const areaScalarField& h,
            const areaScalarField& rho,
            const areaScalarField& p,
            const areaScalarField& Cv
        );


    // Constructors

        //- Construct from components
        frictionModel
        (
            const word& type,
            const dictionary& frictionProperties,
            const areaVectorField& Us,
            const areaScalarField& h,
            const areaScalarField& rho,
            const areaScalarField& p,
            const areaScalarField& Cv
        );


    //- Destructor
    virtual ~frictionModel() = default;


    // Member Functions

        //- Read frictionProperties dictionary
        virtual bool read(const dictionary& frictionProperties) = 0;

        //- Reading field parameters
        void readfield(areaScalarField& fp, word fpName);

        //- Return the friction properties dictionary
        const dictionary& frictionProperties() const
        {
            return frictionProperties_;
        }

        //- Const access to the model coefficients dictionary
        virtual const dictionary& coeffDict() const
        {
            return coeffDict_;
        }

        //- The height limiter
        dimensionedScalar h0() const { return h0_; }

        //- The velocity limiter
        dimensionedScalar u0() const { return u0_; }

        //- Return explicit coefficients for the linear system of equations
        //- considering friction
        virtual const areaVectorField& tauSc() const = 0;

        //- Return implicit coefficients for the linear system of equations
        //- considering friction
        virtual const areaScalarField& tauSp() const = 0;

        //- Return the flow height used as boundary condition for the boundary patch
        virtual const List<scalar>& ht(const float flux, const List<scalar> wide, const List<scalar> slope, const List<scalar> height, const float rhot, const word bdName) const = 0;

        //- Return the uniform flow height used as boundary condition for the boundary patch
        virtual const scalar& htu(const float flux, const List<scalar> wide, const List<scalar> slope, const float rhot, const word bdName) const = 0;

        //- Return the flow velocity used as boundary condition for the boundary patch
        virtual const scalar& ut(const float flux, const List<scalar> wide, const List<scalar> slope, const float rhot, const word bdName) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
