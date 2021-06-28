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

Application
    debrisfaSavageHutterFoam

Description
    Solver for density variable fluids. The equations are depth integrated
    and the finite area method is used.

Author
    Álvaro González Bilbao 

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"
#include "frictionModel.H"
#include "entrainmentModel.H"
#include "depositionModel.H"
#include "SolverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "(debris flows)\n"
        "A depth-integrated solver for debris flows using"
        " Finite Area Methods."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFaMesh.H"
    #include "readGravitationalAcceleration.H"  
    #include "createFaFields.H"
    #include "readTransportProperties.H"
    #include "readFlowCondition.H"
    #include "createTerrainModificationFields.H"
    #include "createTimeControls.H"

    Info<< "\nStarting time loop\n" << endl;

    #include "readSolutionControls.H"

    Info<< nl
        << "Numerical settings" << nl
        << "    max number of iterations " << nCorr << nl
        << "    min number of iterations " << minCorr << nl
        << "    TOL rhoh " << rhohResidualMax << nl
        << "    TOL Cwh " << CwhResidualMax << nl
        << "    TOL Us " << UsResidualMax << nl << endl;

    const bool initDeltaT = runTime.controlDict().get<bool>("initDeltaT");

    if (initDeltaT)
    {
        Info<< "Initializing Delta T" << endl;
        #include "readTimeControls.H"
        #include "surfaceCourantNo.H"
        runTime.setDeltaT
	(
            min(maxCo/(CoNum + SMALL)*runTime.deltaT().value(), maxDeltaT)
        );
    }

    bool final = false;

    while (runTime.run())
    {
        #include "readSolutionControls.H"
        #include "readTimeControls.H"
        #include "surfaceCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

	#include "createFlowCondition.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;
        final = false;

        for (int iCorr = 0; ; iCorr++)
        {
            #include "calcBasalstress.H"

            const areaVectorField & tauSc = friction->tauSc();
            const areaScalarField & tauSp = friction->tauSp();

            faVectorMatrix UsEqn
            (
                fam::ddt(rhoh, Us)
              + shi*fam::div(phi2s, Us)
              + tauSc
              + fam::Sp(tauSp, Us)
             ==
                rho*gs*h
              - fac::grad(pb*h/2.)              
            );
	    
            if (!final)
                UsEqn.relax();

            SolverPerformance<vector> UsResidual = solve(UsEqn);

            tau  = (tauSc + tauSp*Us);
            phis = (fac::interpolate(Us) & aMesh.Le());

            const areaScalarField & Sd = deposition->Sd();
	    const areaScalarField & Sm = entrainment->Sm();

            faScalarMatrix rhohEqn
            (
                fam::ddt(rhoh)
              + fam::div(phis, rhoh)
             ==
                rho_b*Sm
              - fam::Sp
                (
                    rho_b*Sd/(rhoh + dimensionedScalar("small", dimLength*dimDensity, SMALL)),
                    rhoh
                )
            );

            if (!final)
                rhohEqn.relax();

	    solverPerformance rhohResidual = rhohEqn.solve();

            faScalarMatrix CwhEqn
            (
                fam::ddt(Cwh)
              + fam::div(phis,Cwh)
             ==
                Cw_b*Sm
              - fam::Sp
                (
                    Cw_b*Sd/(Cwh + dimensionedScalar("small", dimLength, SMALL)),
                    Cwh
                )
            );

            if (!final)
                CwhEqn.relax();

	    solverPerformance CwhResidual = CwhEqn.solve();

            phi2s = rhohEqn.flux();

            Cw = rho_s/(rhoh/(Cwh+dimensionedScalar("small", dimLength, SMALL))+rho_s-rho_w);
            Cw = min(dimensionedScalar("one", dimless,1),max(dimensionedScalar("zero", dimless,0),Cw));
            h  = max(Cwh/(Cw+dimensionedScalar("small", dimless, SMALL)),hmin);
            rho  = rho_w*Cw+rho_s*(1-Cw);
            Cwh  = Cw*h;
            rhoh = rho*h;
            Cv = 1-Cw;

            Q = (phi2s+CwhEqn.flux()*(rho_s-rho_w))/rho_s;

            deltahh = (Sd-Sm)*runTime.deltaT();

            h.correctBoundaryConditions();
            rhoh.correctBoundaryConditions();
            Cwh.correctBoundaryConditions();
            Us.correctBoundaryConditions();

            if (final)
            {
                if (rhohResidual.initialResidual() < rhohResidualMax && CwhResidual.initialResidual() < CwhResidualMax && mag(UsResidual.initialResidual()) < UsResidualMax)
                {
                    Info<< "reached residual in rhoh = "
                        << rhohResidual.initialResidual()
                        << " < " << rhohResidualMax
			            << ", in Cwh = "
                        << CwhResidual.initialResidual()
                        << " < " << CwhResidualMax
                        << " and in Us = "
                        << UsResidual.initialResidual()
                        << " < " << UsResidualMax
                        << ", stopping loop!" << endl;
                }
                else
                {
                    Info<< "Reached maximum numbers of iterations, "
                        << "stopping loop!" << endl;
                }
                break;
            }

            if
            (
                (
                    rhohResidual.initialResidual() < rhohResidualMax
                 && CwhResidual.initialResidual() < CwhResidualMax
                 && mag(UsResidual.initialResidual()) < UsResidualMax
                 && iCorr >= minCorr
                )
             || iCorr >= nCorr
            )
            {
                final = true;
            }
        }

	he += deltahh;

        if (terrainModification)
	{  
		deltah += deltahh;

		#include "terrainModification.H"

		deltac0 = (c & vector(0,0,1))-c0;
	}
	else
	{
		deltah0 = he-he0;
	}

        if (runTime.outputTime())
        {
            runTime.write();

            if (terrainModification)
            {
                deltac0.write();	
            }
            else
            {
                deltah0.write();
            }
        }

        runTime.printExecutionTime(Info);
    }
   
    Info<< nl << "End" << endl;

    return 0;
}

// ************************************************************************* //
