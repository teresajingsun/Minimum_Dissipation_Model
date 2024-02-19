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

#include "scalarQR.H"		
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam		
{			
namespace LESModels	
{			

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>		
tmp<volScalarField> scalarQR<BasicTurbulenceModel>::k
(			
    const tmp<volTensorField>& gradU	
) const
{                                             
    volSymmTensorField  D(symm(gradU()));       
    const volVectorField& U = this->U_;
    const auto& rhok = U.mesh().lookupObject<volScalarField>("rhok");
    volScalarField	buoyancy(1.0/4.0 * I & g_ & gradU & fvc::grad(rhok));

    volScalarField 	r(-det(D));
    volScalarField      q(max((1.0/2.0)*dev(D) && D, dimensionedScalar("q",dimensionSet(0, 0, -2, 0, 0, 0, 0), 1e-15)));	
    volScalarField	rBuoyancy(max(r + buoyancy, dimensionedScalar("buoyancy", dimensionSet(0, 0, -3, 0, 0, 0, 0), 0.0))); 
    
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
		sqr(this->delta()*mag(rBuoyancy)/q)								
	)
    );	
}   //end const           	
                	

template<class BasicTurbulenceModel>
void scalarQR<BasicTurbulenceModel>::correctNut()
{               	
    volScalarField k(this->k(fvc::grad(this->U_)));

    this->nut_ = Ck_*this->delta()*sqrt(k);		//correctNut is eddy viscosity. scalarQR's Ck is equal to Ck^2 in Smagorinsky.C
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*
The template parameter of scalarQR class is BasicTurbulenceModel. 
The template parameter of LESeddyViscosity class is BasicTurbulenceModel.
scalarQR inherits from LESeddyViscosity base class.
The parameterized constructor of the base class can only be called using initializer list, here is an exampel
    B::B(type x):A(x){}   where the derived class is B, the base class is A.   
*/

template<class BasicTurbulenceModel>
scalarQR<BasicTurbulenceModel>::scalarQR			//Constructor function has no return values neither constructor prototype declaration within the class nor constructor definition.
(						//Constructor simply initialize the object
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
    LESeddyViscosity<BasicTurbulenceModel>      //calling the parameterized constructor of Base class LESeddyViscosisty using Initializer list
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

    //QR's Ck is comparable to Ck^2 in Smagorinsky.C
    Ck_                                       
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            0.024
        )
    ),
    g_
    (
        "g",
        dimLength/sqr(dimTime),
        meshObjects::gravity::New(this->mesh_.time()).value()
    )

/*
    gh_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gh",
            this->coeffDict_,
            dimLength/sqr(dimTime),
            9.81
        )
    )

    beta_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta",
            this->coeffDict_,
            dimless/dimTemperature,
            3.3e-03
        )
    ),

    TRef_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "TRef",
            this->coeffDict_,
            dimTemperature,
	    300
        )
    ),
    g_(meshObjects::gravity::New(this->mesh_.time()))

*/

{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool scalarQR<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        Ck_.readIfPresent(this->coeffDict());
	//beta_.readIfPresent(this->coeffDict());
	//TRef_.readIfPresent(this->coeffDict());
	//gh_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void scalarQR<BasicTurbulenceModel>::correct()
{
    LESeddyViscosity<BasicTurbulenceModel>::correct();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
