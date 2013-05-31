//Accessing nu1 and nu2
dimensionedScalar nu1(transportProperties.subDict("phase1").lookup("nu"));
dimensionedScalar nu2(transportProperties.subDict("phase2").lookup("nu"));

//Accessing rho1 and rho2
dimensionedScalar rho1f(transportProperties.subDict("phase1").lookup("rho"));
dimensionedScalar rho2f(transportProperties.subDict("phase1").lookup("rho"));

//Converting volScalarField to surfaceScalarField
surfaceScalarField alpha1f(fvc::interpolate(alpha1));


//Arithmetic Mean
surfaceScalarField muf1 = alpha1f*rho1f*nu1 + (scalar(1) - alpha1f)*rho2f*nu2;

//Harmonic Mean
surfaceScalarField muf2 = rho1f*nu1*rho2f*nu2/(rho1f*nu1*(scalar(1) - alpha1f) + rho2f*nu2*alpha1f);

//Unit interface Normal Flux 
surfaceScalarField nf = nHatfv & mesh.Sf();

//Interpolated KotheViscosity
surfaceScalarField muf=nf*muf1+(scalar(1)-nf)*muf2;
