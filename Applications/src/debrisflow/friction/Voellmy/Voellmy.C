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

#include "Voellmy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace frictionModels
{
    defineTypeNameAndDebug(Voellmy, 0);
    addToRunTimeSelectionTable(frictionModel, Voellmy, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::frictionModels::Voellmy::Voellmy
(
    const dictionary& frictionProperties,
    const areaVectorField& Us,
    const areaScalarField& h,
    const areaScalarField& rho,
    const areaScalarField& p,
    const areaScalarField& Cv
)
:
    frictionModel(type(), frictionProperties, Us, h, rho, p, Cv),
    mu_
    (
        IOobject
        (
            "mu",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimless)
    ),
    xi_
    (
        IOobject
        (
            "xi",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimAcceleration)
    )
{
    readfield(mu_,"mu");
    readfield(xi_,"xi");
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


const Foam::areaVectorField& Foam::frictionModels::Voellmy::tauSc() const
{
    resetTauSc();

    return tauSc_;
}


const Foam::areaScalarField& Foam::frictionModels::Voellmy::tauSp() const
{
    resetTauSp();

    areaScalarField u(mag(Us_));
    areaScalarField lambda
    (
      rho_*dimensionedScalar("g", dimAcceleration, 9.81)/(xi_+dimensionedScalar("small", dimAcceleration, SMALL))
    );

    // Coulomb Part
    tauSp_ += p_ * mu_ *
              1./(u + u0_);

    // Kinetic Part
    tauSp_ += lambda * u;

    return tauSp_;
}


bool Foam::frictionModels::Voellmy::read
(
    const dictionary& frictionProperties
)
{
    readDict(type(), frictionProperties);

    return true;
}

const Foam::scalar& Foam::frictionModels::Voellmy::ut
(
    const float flux, 
    const List<scalar> wide,
    const List<scalar> slope,
    const float rhot,
    const word bdName
) const
{
    if (wide.size() != 0 && flux != 0)
    {
	    faBoundaryMesh faBoundaries
	    (
		IOobject
		(
		    "faBoundary",
		    Us_.mesh().time().constant(),
		    "faMesh",
		    Us_.mesh().mesh(),
		    IOobject::MUST_READ,
		    IOobject::NO_WRITE
		),
		Us_.mesh()
	    );

	    label index = faBoundaries.findIndex(bdName);
	    const faPatch faPatchDict(faBoundaries[index],faBoundaries);
	    const labelUList edgeFaces = faPatchDict.edgeFaces();

	    float mu = 0.;
	    float xi = 0.;
	
	    forAll(edgeFaces, edgeI)
	    {
		mu += mu_[edgeFaces[edgeI]];
		xi += xi_[edgeFaces[edgeI]];
	    };

	    mu /= edgeFaces.size();
	    xi /= edgeFaces.size();

	    float lambda;
	    float sum = 0.;
	    float B = 0.;

	    forAll(wide, I)
	    {
		B += wide[I];
		sum += wide[I]*pow(max(sin(slope[I]*M_PI/180)-mu*cos(slope[I]*M_PI/180),0),1/2.);
	    };

	    lambda = pow(sum, 2/3.);
	    ut_ = pow(flux*xi,1/3.)/B*lambda;
    }
    else
    {
	    ut_ = 0;
    }

    return ut_;
}

const Foam::scalar& Foam::frictionModels::Voellmy::htu
(
    const float flux, 
    const List<scalar> wide,
    const List<scalar> slope,
    const float rhot,
    const word bdName
) const
{
    scalar htu;

    if (wide.size() != 0 && flux != 0)
    {
	    faBoundaryMesh faBoundaries
	    (
		IOobject
		(
		    "faBoundary",
		    Us_.mesh().time().constant(),
		    "faMesh",
		    Us_.mesh().mesh(),
		    IOobject::MUST_READ,
		    IOobject::NO_WRITE
		),
		Us_.mesh()
	    );

	    label index = faBoundaries.findIndex(bdName);
	    const faPatch faPatchDict(faBoundaries[index],faBoundaries);
	    const labelUList edgeFaces = faPatchDict.edgeFaces();

	    float mu = 0.;
	    float xi = 0.;

	    forAll(edgeFaces, edgeI)
	    {
		mu += mu_[edgeFaces[edgeI]];
		xi += xi_[edgeFaces[edgeI]];
	    };

	    mu /= edgeFaces.size();
	    xi /= edgeFaces.size();

	    float lambda;
	    float sum = 0.;

	    forAll(wide, I)
	    {
		sum += wide[I]*pow(max(sin(slope[I]*M_PI/180)-mu*cos(slope[I]*M_PI/180),0),1/2.);
	    };

	    lambda = pow(sum, 2/3.);
	    htu = pow(flux/pow(xi, 1/2.),2/3.)/lambda;
    }
    else
    {
	    htu = 0.;
    }

    return htu;
}

const Foam::List<Foam::scalar>& Foam::frictionModels::Voellmy::ht
(
    const float flux, 
    const List<scalar> wide,
    const List<scalar> slope,
    const List<scalar> height,
    const float rhot,
    const word bdName
) const
{
    if (wide.size() != 0 && flux != 0)
    {
	    List<scalar> htt(height.size(), zero());

	    float sum_wide = 0.;
	    float sum_height = 0.;
	    float z_top;
	    scalar utt = ut(flux, wide, slope, rhot, bdName);

	    forAll(wide, I)
	    {
		sum_wide += wide[I];
		sum_height += height[I]*wide[I];
	    };

	    z_top = ((flux/utt)+sum_height)/sum_wide;

	    forAll(htt, I)
	    {
		htt[I] = max(z_top-height[I], 0.);
	    };

	    bool final = true;
	    float new_sum_wide = 0.;
	    int nIter = 0;

	    while (final == true && nIter<5)
	    {
		    new_sum_wide = sum_wide;
		    sum_wide = 0.;
		    sum_height = 0.;

		    forAll(wide, I)
		    {
			if (htt[I]>0)
			{
				sum_wide += wide[I];
				sum_height += height[I]*wide[I];
			}
		    };

		    z_top = ((flux/utt)+sum_height)/sum_wide;

		    forAll(htt, I)
		    {
			htt[I] = max(z_top-height[I], 0.);
		    };

		    if (new_sum_wide == sum_wide)
		    {
			final = false;
		    }
		    else
		    {
			nIter ++;
		    }
	    }

	    if (nIter == 5)
	    {
		    scalar htuu = htu(flux, wide, slope, rhot, bdName);
		    List<scalar> htt(height.size(),  htuu);
	    }

	    ht_ = htt;
    }
    else
    {
	    List<scalar> htt(height.size(), zero());
	    ht_ = htt;
    }

    return ht_;
}

// ************************************************************************* //
