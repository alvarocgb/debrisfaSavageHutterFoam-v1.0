
    scalar nCorr =
            aMesh.solutionDict().lookupOrDefault<int>("nOuterCorrectors", 50);

    scalar minCorr =
            aMesh.solutionDict().lookupOrDefault<int>("minCorrectors", 3);

    scalar rhohResidualMax =
        aMesh.solutionDict().lookupOrDefault<scalar>("rhohResidualMax", 0.);

    scalar CwhResidualMax =
        aMesh.solutionDict().lookupOrDefault<scalar>("CwhResidualMax", 0.);

    scalar UsResidualMax =
        aMesh.solutionDict().lookupOrDefault<scalar>("UsResidualMax", 1);

