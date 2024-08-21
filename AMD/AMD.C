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

#include "AMD.H"
#include "fvOptions.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class BasicTurbulenceModel>
void AMD<BasicTurbulenceModel>::calcFilter()
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

	    filter_[celli].x() = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
	    filter_[celli].y() = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
	    filter_[celli].z() = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));
	    // And similar for yDim and zDim
	    //Info << filter_[celli].x() << " " << filter_[celli].y() << " " << filter_[celli].z() << endl;
	  }
	}	

template<class BasicTurbulenceModel>
tmp<volScalarField> AMD<BasicTurbulenceModel>::k
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
void AMD<BasicTurbulenceModel>::correctNut()
{
    volTensorField gradU = fvc::grad(this->U_);
    volScalarField xDim=filter_.component(vector::X);
    volScalarField yDim=filter_.component(vector::Y);
    volScalarField zDim=filter_.component(vector::Z);
    //First row
    volScalarField gradUxx = gradU.component(tensor::XX);
    volScalarField gradUxy = gradU.component(tensor::XY);
    volScalarField gradUxz = gradU.component(tensor::XZ);

    //Second row
    volScalarField gradUyx = gradU.component(tensor::YX);
    volScalarField gradUyy = gradU.component(tensor::YY);
    volScalarField gradUyz = gradU.component(tensor::YZ);

    //Third row
    volScalarField gradUzx = gradU.component(tensor::ZX);
    volScalarField gradUzy = gradU.component(tensor::ZY);
    volScalarField gradUzz = gradU.component(tensor::ZZ);

    volScalarField numerator_xx = sqr(xDim) * gradUxx * ( gradUxx*gradUxx + gradUyx*gradUxy + gradUzx*gradUxz);
    volScalarField numerator_xy = sqr(xDim) * gradUyx * ( gradUxx*gradUyx + gradUyx*gradUyy + gradUzx*gradUyz);
    volScalarField numerator_xz = sqr(xDim) * gradUzx * ( gradUxx*gradUzx + gradUyz*gradUzy + gradUzx*gradUzz);

    volScalarField numerator_yx = sqr(yDim) * gradUxy * ( gradUxy*gradUxx + gradUyy*gradUxy + gradUzy*gradUxz);
    volScalarField numerator_yy = sqr(yDim) * gradUyy * ( gradUxy*gradUyx + gradUyy*gradUyy + gradUzy*gradUyz);
    volScalarField numerator_yz = sqr(yDim) * gradUzy * ( gradUxy*gradUzx + gradUyy*gradUzy + gradUzy*gradUzz);

    volScalarField numerator_zx = sqr(zDim) * gradUxz * ( gradUxz*gradUxx + gradUyz*gradUxy + gradUzz*gradUxz);
    volScalarField numerator_zy = sqr(zDim) * gradUyz * ( gradUxz*gradUyx + gradUyz*gradUyy + gradUzz*gradUyz);
    volScalarField numerator_zz = sqr(zDim) * gradUzz * ( gradUxz*gradUzx + gradUyz*gradUzy + gradUzz*gradUzz);

//    Info << numerator_xx << " " << numerator_xy << " " << numerator_xz 
//	 << numerator_yx << " " << numerator_yy << " " << numerator_yz 
//    	 << numerator_zx << " " << numerator_zy << " " << numerator_zz 
//	 << endl;

    volScalarField numerator
	    	    (max 
		     	( -  (numerator_xx + numerator_xy + numerator_xz +
			      numerator_yx + numerator_yy + numerator_yz +
			      numerator_zx + numerator_zy + numerator_zz),
		    	  dimensionedScalar("numerator", dimensionSet(0, 2, -3, 0, 0, 0, 0), 1e-25)
			)
		    );
    
    this->nut_ = Ck_ * numerator / (gradU && gradU);
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
AMD<BasicTurbulenceModel>::AMD
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
	  dimensionedVector ("filter",dimLength,vector::zero)
     )

{
    //Info << "AMD Constructor" << endl;
    if (type == typeName)
    {
	//Info << "after typename check" << endl;
	calcFilter();
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool AMD<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        Ck_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void AMD<BasicTurbulenceModel>::correct()
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
