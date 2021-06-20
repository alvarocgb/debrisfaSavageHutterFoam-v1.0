pb -= shi*fac::ndiv(phi2s, Us)&n;
pb += rho*(g*h)&n;

if (pressureFeedback)
{
    pb -= fac::ngrad(pb*h/2.)&n;
}

pb = Foam::max(pb, dimensionedScalar(dimPressure));