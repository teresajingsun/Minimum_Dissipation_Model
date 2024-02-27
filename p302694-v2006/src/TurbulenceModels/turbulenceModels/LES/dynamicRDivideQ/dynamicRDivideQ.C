/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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



dynamicRDivideQ - Implementation of the dynamic Smagorinsky
                     SGS model.
Copyright Information
    Copyright (C) 1991-2009 OpenCFD Ltd.
    Copyright (C) 2010-2021 Alberto Passalacqua 
License
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "dynamicRDivideQ.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void dynamicRDivideQ<BasicTurbulenceModel>::correctNut
(
    const tmp<volTensorField>& gradU
)
{
    const volSymmTensorField S(symm(gradU));
    const volScalarField     q(max((1.0/2.0)*(dev(S) && S), dimensionedScalar("q", dimensionSet(0, 0, -2, 0, 0, 0, 0), 1e-10)));
    const volScalarField     r(max(-det(S), dimensionedScalar("r",dimensionSet(0, 0, -3, 0, 0, 0, 0), 0.0)));
    
    // The SGS viscosity is bounded so that eddy viscosity cannot be equal to or less than zero.
    // No warning message is printed when this limitation is applied.
  
    this->nut_ = max(cD_*sqr(this->delta())*mag(r)/q, dimensionedScalar("nut", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0));  
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
void dynamicRDivideQ<BasicTurbulenceModel>::correctNut()
{
    correctNut(fvc::grad(this->U_));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
dynamicRDivideQ<BasicTurbulenceModel>::dynamicRDivideQ
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    LESeddyViscosity<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    cD_
    (
        IOobject
        (
            IOobject::groupName("cD_", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("cD", dimless, 0.0)
    ),

    cI_
    (
        IOobject
        (
            IOobject::groupName("cI_", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("cI", dimless, 0.0)
    ),
    filterPtr_(LESfilter::New(U.mesh(), this->coeffDict())),
    filter_(filterPtr_())
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void dynamicRDivideQ<BasicTurbulenceModel>::calcCD
(
    const volSymmTensorField& S
)
{
    const volVectorField& U = this->U_;
    tmp<volSymmTensorField> LL = dev(filter_(sqr(U)) - (sqr(filter_(U))));
    const volScalarField  q
    (
    	max
	(
	  dimensionedScalar("q1", dimensionSet(0, 0, -2, 0, 0, 0, 0), 1e-10)
	  , (1.0/2.0)*(dev(S) && S)
	)
    );
    const volSymmTensorField MM
    (
        sqr(this->delta())
       *(filter_((mag(-det(S))/q)*(S)) 
       - 4.0*(mag((-det(filter_(S))))/((1.0/2.0)*(filter_(dev(S)) && filter_(S))))*filter_(dev(S)))
    );

    // Locally averaging MMMM on cell faces
    volScalarField MMMM = fvc::average(magSqr(MM));

    MMMM.max(VSMALL);

    // Performing local average on cell faces on return
    cD_ = -0.5*fvc::average(LL && MM)/MMMM;

    if (debug)
    {
        Info<< "min(cD) = " << min(cD_) 
            << ", max(cD) = " << max(cD_)
            << ", average(cD) = " << average(cD_) << endl;
    }
}

template<class BasicTurbulenceModel>
void dynamicRDivideQ<BasicTurbulenceModel>::calcCI
(
    const volSymmTensorField& S
)
{
    const volVectorField& U = this->U_;
    tmp<volScalarField> KK = 0.5*(filter_(magSqr(U)) - magSqr(filter_(U)));

    const volScalarField mm
    (
        sqr(this->delta())*(4*sqr(mag(filter_(dev(S)))) - filter_(sqr(mag(dev(S)))))
    );

    // Locally averaging mmmm on cell faces
    volScalarField mmmm = fvc::average(magSqr(mm));

    mmmm.max(VSMALL);

    // Performing local average on cell faces on return
    cI_ = fvc::average(KK*mm)/mmmm;

    if (debug)
    {
        Info<< "min(cI) = " << min(cI_) 
            << ", max(cI) = " << max(cI_)
            << ", average(cI) = " << average(cI_) << endl;
    }
}


template<class BasicTurbulenceModel>
bool dynamicRDivideQ<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        filter_.read(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void dynamicRDivideQ<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField S(dev(symm(gradU)));

    calcCD(S);
    calcCI(S);
    correctNut(gradU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
