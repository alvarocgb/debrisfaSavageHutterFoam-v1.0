IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

autoPtr<frictionModel> friction
(
    frictionModel::New(transportProperties, Us, h, rho, pb, Cv)
);

autoPtr<entrainmentModel> entrainment
(
    entrainmentModel::New(transportProperties, Us, h, he, rho, pb, tau, Cv)
);

autoPtr<depositionModel> deposition
(
    depositionModel::New(transportProperties, Us, h, he, rho, pb, tau, Cv)
);

dimensionedScalar hmin("hmin", dimLength, transportProperties);
dimensionedScalar shi("shi", dimless, transportProperties);
dimensionedScalar rho_w("rho_w", dimDensity, transportProperties);
dimensionedScalar rho_s("rho_s", dimDensity, transportProperties);
dimensionedScalar rho_b("rho_b", dimDensity, transportProperties);
dimensionedScalar Cw_b("Cw_b", dimless, (rho_s.value()-rho_b.value())/(rho_s.value()-rho_w.value()));

Switch pressureFeedback(transportProperties.lookup("pressureFeedback"));
Switch terrainModification(transportProperties.lookup("terrainModification"));

Info << "Running with" << endl
     << "    hmin " << hmin << endl
     << "    rho_w " << rho_w << endl
     << "    rho_s " << rho_s << endl
     << "    rho_b " << rho_b << endl
     << "    Cw_b " << Cw_b << endl
     << "    pressureFeedback is " << pressureFeedback << endl
     << "    terrainModification is " << terrainModification << endl;

int nTerrainIter = 0;
dimensionedScalar deltahMin("deltahMin", dimLength, 0);
dimensionedScalar hMinterrain("hMinterrain", dimLength, 0);
scalar alpha_relax = 0;
scalar alpha_smooth = 0;
scalar alpha_gamma = 0;
scalar beta_gamma = 0;

if (terrainModification)
{
     dictionary terrainModificationCoeffs = transportProperties.optionalSubDict("terrainModificationCoeffs");
     nTerrainIter = terrainModificationCoeffs.lookupOrDefault<int>("nTerrainIter",5);
     deltahMin = dimensionedScalar("deltahMin", dimLength, terrainModificationCoeffs);
     hMinterrain = dimensionedScalar("hMinterrain", dimLength, terrainModificationCoeffs);
     alpha_relax = terrainModificationCoeffs.lookupOrDefault<scalar>("alpha_relax",0.85);
     alpha_smooth = terrainModificationCoeffs.lookupOrDefault<scalar>("alpha_smooth",0.7);
     alpha_gamma = terrainModificationCoeffs.lookupOrDefault<scalar>("alpha_gamma",0.1);
     beta_gamma = terrainModificationCoeffs.lookupOrDefault<scalar>("beta_gamma",0.6);
     Info << "         deltahMin " << deltahMin << endl
          << "         hMinterrain " << hMinterrain << endl
          << "         nTerrainIter " << nTerrainIter << endl
          << "         alpha_relax " << alpha_relax << endl
          << "         alpha_smooth " << alpha_smooth << endl
          << "         alpha_gamma " << alpha_gamma << endl 
          << "         beta_gamma " << beta_gamma << nl << endl; 
}
