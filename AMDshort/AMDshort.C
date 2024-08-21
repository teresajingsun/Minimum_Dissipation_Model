/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "AMDshort.H"
#include "fvOptions.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicTurbulenceModel>
void AMDshort<BasicTurbulenceModel>::calcFilter()
	{
	  //Info << "calculate filter" <<endl;
	  const faceList & ff = this->mesh_.faces();

	  const pointField & pp = this->mesh_.points();

	  forAll ( this->mesh_.C(), celli)
	  {
	    const cell & cc = this->mesh_.cells()[celli];

	    labelList pLabels(cc.labels(ff));

	    pointField pLocal(pLabels.size(), vector::zero);

	    forAll (pLabels, pointi)
	        pLocal[pointi] = pp[pLabels[pointi]];

	    vector tmp;
	    tmp.x() = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
	    tmp.y() = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
	    tmp.z() = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));
	    
	    filter_[celli].diag(tmp);
	    //Info << filter_[celli].x() << " " << filter_[celli].y() << " " << filter_[celli].z() << endl;
	  }
	}	

template<class BasicTurbulenceModel>
tmp<volScalarField> AMDshort<BasicTurbulenceModel>::k
(
    const tmp<volTensorField>& gradU
) const
{
    volSymmTensorField D(symm(gradU));
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("k", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
	    
	    (2.0*Ck_/this->Ce_)*sqr(this->delta())*magSqr(dev(D))
        )
    );
}


template<class BasicTurbulenceModel>
void AMDshort<BasicTurbulenceModel>::correctNut()
{

    volTensorField gradU = fvc::grad(this->U_);

    volTensorField numerator1 ( filter_ & gradU ); 
    
    volTensorField numerator2 ( (filter_ & gradU) & (filter_ & gradU) ); 
    
    volScalarField numerator3 ( gradU.T() & (filter_.T() ) & (filter_ & gradU) && symm(gradU) ); 

    volScalarField numerator  ( max ( - numerator3,
		    	  dimensionedScalar("numerator", dimensionSet(0, 2, -3, 0, 0, 0, 0), 1e-25) ) );

    volScalarField denominator (gradU && gradU);
    
    this->nut_ = Ck_* numerator/ denominator;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
AMDshort<BasicTurbulenceModel>::AMDshort
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

    Ck_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            0.094
        )
    ),

    filter_
    (
     	IOobject
	  (
	      "filter",
	      U.mesh().time().timeName(),
	      U.mesh(),
	      IOobject::NO_READ,
	      IOobject::AUTO_WRITE
	  ),
	  U.mesh(),
	  dimensionedTensor ("filter",dimLength,tensor::zero)
     )

{
    //Info << "AMDshort Constructor" << endl;
    if (type == typeName)
    {
	//Info << "after typename check" << endl;
	calcFilter();
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool AMDshort<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        Ck_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void AMDshort<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    LESeddyViscosity<BasicTurbulenceModel>::correct();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
