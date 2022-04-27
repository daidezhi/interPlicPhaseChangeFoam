{
    // Temporarily making U relative to mesh motion
    if (mesh.moving())
    {
        U -= fvc::reconstruct(mesh.phi());
    }

    plicVofSolver.preProcess();

    plicVofSolver.calcOrientation();
    plicVofSolver.reconstruct();

    zeroField Sp;
    zeroField Su;
    plicVofSolver.advect(Sp, (Su + (rho2/(rho2-rho1))*divU)());
    //plicVofSolver.advect((Sp + vDotvmcAlphal)(), (Su + vDotcAlphal)());


    // Making U absolute again after advection step
    if (mesh.moving())
    {
        U += fvc::reconstruct(mesh.phi());
    }

    rhoPhi = plicVofSolver.getRhoPhi(rho1, rho2);

    //Info<< "plicVofSolving: Mass conservation Error = "
    //    << plicVofSolver.massConservationError() << endl;

    alpha2 = 1.0 - alpha1;
    mixture->correct();
}

Info<< "Phase-1 volume fraction = "
    << alpha1.weightedAverage(mesh.Vsc()).value()
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
    << endl;
