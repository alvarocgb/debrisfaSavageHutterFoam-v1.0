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
    Álvaro González Bilbao

\*---------------------------------------------------------------------------*/

#include "QuadraticManning.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace frictionModels
{
    defineTypeNameAndDebug(QuadraticManning, 0);
    addToRunTimeSelectionTable(frictionModel, QuadraticManning, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::frictionModels::QuadraticManning::QuadraticManning
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
    g_("g", dimAcceleration, 9.81),
    K_
    (
        IOobject
        (
            "K",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimless)
    ),
    n_
    (
        IOobject
        (
            "n",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimless)
    ),
    alpha1_
    (
        IOobject
        (
            "alpha1",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimDynamicViscosity)
    ),
    alpha2_
    (
        IOobject
        (
            "alpha2",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimPressure)
    ),
    beta1_
    (
        IOobject
        (
            "beta1",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimless)
    ),
    beta2_
    (
        IOobject
        (
            "beta2",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimless)
    ),
    Ty_
    (
        IOobject
        (
            "Ty",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimPressure)
    ),
    eta_
    (
        IOobject
        (
            "eta",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimDynamicViscosity)
    ),
    nM_
    (
        IOobject
        (
            "n",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimless)
    ),
    Cvlim_
    (
        IOobject
        (
            "Cvlim",
            Us_.time().timeName(),
            Us_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Us_.mesh(),
        dimensionedScalar(dimless)
    )
{
    readfield(K_,"K");
    readfield(n_,"n");
    readfield(alpha1_,"alpha1");
    readfield(alpha2_,"alpha2");
    readfield(beta1_,"beta1");
    readfield(beta2_,"beta2");
    readfield(nM_,"nM");
    readfield(Cvlim_,"Cvlim");
}

// * * * * * * * * * * * * * Protected Functions  * * * * * * * * * * * * * * //

void Foam::frictionModels::QuadraticManning::resetTy() const
{
    Ty_ = dimensionedScalar(dimPressure);
}

void Foam::frictionModels::QuadraticManning::reseteta() const
{
    eta_ = dimensionedScalar(dimDynamicViscosity);
}

void Foam::frictionModels::QuadraticManning::calcTy() const
{
    resetTy();

    Ty_ = alpha2_*exp(Cv_*beta2_);
}

void Foam::frictionModels::QuadraticManning::calceta() const
{
    reseteta();

    eta_ = alpha1_*exp(Cv_*beta1_);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


const Foam::areaVectorField& Foam::frictionModels::QuadraticManning::tauSc() const
{
    resetTauSc();

    return tauSc_;
}


const Foam::areaScalarField& Foam::frictionModels::QuadraticManning::tauSp() const
{
    resetTauSp();
    calcTy();
    calceta();

    areaScalarField u(mag(Us_));
    areaScalarField alpha = pos(Cv_-Cvlim_);

    // Yield Stress
    tauSp_ += Ty_ *
              1./(u + u0_)*alpha;

    // Viscous stress
    tauSp_ += 1./(8*(h_+h0_)) * K_ * eta_*alpha;

    // Turbulent-dispersive stress
    tauSp_ += rho_*n_*n_*g_.value()*u/pow((h_+h0_)*dimensionedScalar("n2", dimAntiL, 1.0), 1/3.)*alpha;

    tauSp_ += (rho_*nM_*nM_*g_.value()*(u+u0_)/pow((h_+h0_)*dimensionedScalar("n2", dimAntiL, 1.0), 1/3.))*(dimensionedScalar("one", dimless,1)-alpha);

    return tauSp_;
}

bool Foam::frictionModels::QuadraticManning::read
(
    const dictionary& frictionProperties
)
{
    readDict(type(), frictionProperties);

    return true;
}

const Foam::scalar& Foam::frictionModels::QuadraticManning::ut
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
	    if (flux == 0)
	    {
		ut_ = 0;
	    }
	    else
	    {
		scalar htuu = htu(flux, wide, slope, rhot, bdName);

		float B = 0.;

		forAll(wide, I)
		{	
			B += wide[I];
		};  

		ut_ = flux/(B*htuu);
	    }
    }
    else
    {
	    ut_ = 0;
    }

    return ut_;
}

const Foam::scalar& Foam::frictionModels::QuadraticManning::htu
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
	    calcTy();
	    calceta();

	    float B = 0.;
	    float avrg_slope = 0.;

	    forAll(wide, I)
	    {
		B += wide[I];
		avrg_slope += slope[I];
	    };

	    avrg_slope /= slope.size();

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

	    float n = 0.;
	    float K = 0.;
	    float eta = 0.;
	
	    forAll(edgeFaces, edgeI)
	    {
		n  += n_[edgeFaces[edgeI]];
		K += K_[edgeFaces[edgeI]];
		eta  += eta_[edgeFaces[edgeI]];
	    };

	    n /= edgeFaces.size();
	    K /= edgeFaces.size();
	    eta /= edgeFaces.size();

	    htu = pow(flux*n/(B*pow(Foam::sin(M_PI/180*avrg_slope),1/2.)),3/5.);
	    float utu = flux/(B*htu);
	    float new_utu;

	    bool final = false; 
	   
	    while (!final)
	    {
		    new_utu = pow(flux/B*g_.value()*Foam::sin(M_PI/180*avrg_slope)/(n*n*g_.value()/pow(htu, 1/3.)+K*eta*B/(8*flux*rhot)),1/3.);

		    if (abs(new_utu-utu)/(SMALL+utu) < 0.01)
		    {
			final = true;
		    }

		    utu = (utu+new_utu)/2;
		    htu = flux/(B*utu);
	    }
    }
    else
    {
	    htu = 0;
    }

    return htu;
}

const Foam::List<Foam::scalar>& Foam::frictionModels::QuadraticManning::ht
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
