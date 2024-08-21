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

\*---------------------------------------------------------------------------*/

#include "QRgrad.H"		
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam		
{			
namespace LESModels	
{			

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


//QR subgrid model do not provide information about the turbulent kinetic energy.
template<class BasicTurbulenceModel>		
tmp<volScalarField> QRgrad<BasicTurbulenceModel>::k 
(			
    const tmp<volTensorField>& gradU	
) const
{                                             
    volSymmTensorField  D(symm(gradU()));       
    volTensorField	product1 (- gradU() & gradU());
    volTensorField	product2 (- gradU().T() & gradU());

    volScalarField	numerator(max( - gradU().T() & gradU() && D, 
			    dimensionedScalar("numerator",dimensionSet(0, 0, -3, 0, 0, 0, 0), 1e-25)));

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
		sqr(this->delta()*mag(numerator) / (gradU && gradU))
	)
      );	
}   //end const           	
                	

template<class BasicTurbulenceModel>
void QRgrad<BasicTurbulenceModel>::correctNut()
{               	
    volScalarField k(this->k(fvc::grad(this->U_)));

    this->nut_ = Ck_*this->delta()*sqrt(k);		//correctNut is eddy viscosity. QRgrad's Ck is equal to Ck^2 in Smagorinsky.C
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class BasicTurbulenceModel>
QRgrad<BasicTurbulenceModel>::QRgrad			
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
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool QRgrad<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        Ck_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void QRgrad<BasicTurbulenceModel>::correct()
{
    LESeddyViscosity<BasicTurbulenceModel>::correct();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
