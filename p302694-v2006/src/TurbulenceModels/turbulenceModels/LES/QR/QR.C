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

#include "QR.H"		
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
namenspace groups named entities that would otherwise have a global scope into narrow scopes. 
LESModels is a normal variable declared within a namespace Foam. Namespce can extend across different files of source code.
To access the variable namespace LESModels from outside Foam they should be qualified like: Foam::namespace LESModels
*/
namespace Foam		
{			
namespace LESModels	
{			

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

/*
template <class identifier> function_declaration | template <typename identifier> function_declaration
BasicTurbuleceModel is defined as a generic type, a template parameter
Calling function template: function_name <type> (function arguements)
*/
template<class BasicTurbulenceModel>		
tmp<volScalarField> QR<BasicTurbulenceModel>::k 
(			
    const tmp<volTensorField>& gradU	
) const
{                                             
    volSymmTensorField  D(symm(gradU));       
/*
Dimension of velocity M^0LT^(-1), dimension of gradient M^0LT^0
Dimension of velocity gradient, i.e.,velocity/distance, equal to M^0L^0T(-1)
In OpenFOAM physical dimensions with respect to the SI unit system: SI [kg, m, s, K, mol, A, cd]
q is the second invariant of velocity gradient. Dimension of q is (velocity gradient)^2, i.e., M^0L^0T^(-2)    
r is the third invariant.  Dimension of r (velocity gradient)^3, i.e., M^0L^0T^(-3)
*/
    volScalarField 	r(max(-det(D),dimensionedScalar("r",dimensionSet(0, 0, -3, 0, 0, 0, 0), 0.0)));
    volScalarField      q(max((1.0/2.0)*dev(D) && D, dimensionedScalar("q",dimensionSet(0, 0, -2, 0, 0, 0, 0), 1e-15)));	
    //Setting the dimension of q1 as same as q, to compare q and q1 in If statement. But can't solve the case 0 < q << 1 that eddy viscosity tends to infinite. 
//    volScalarField	q1(0.0*dev(D) && D);
//
//    if ( q <= q1 )              //lhs and rhs should be the same type
//    {
//    return tmp<volScalarField>
//	(
//	new	volScalarField		
//	(
//		IOobject
//		(
//			IOobject::groupName("k",this->alphaRhoPhi_.group()),
//			this->runTime_.timeName(),
//			this->mesh_
//		),
//		q1		// k=0 when q<=0; Potential issue would be that the dimensions of q1 and k are not match                              
//	)
//      );
// 
//    }	
     return tmp<volScalarField>
      (
	new volScalarField	
	(
		IOobject
		(
			IOobject::groupName("k", this->alphaRhoPhi_.group()),//k is the name of file containing the dictionary
			this->runTime_.timeName(),		//tells OpenFOAM to save the file in a dictionary called as the current time	
			this->mesh_				//pointer_name -> variable_name(a sturcture or union)
								//mesh is the objectRegistry
		),						//Data member of class, both static and non-static, are named like ordinary variable,
		//mag(r)/q					//but with a trailing underscore, like: mesh_					
		sqr(this->delta()*mag(r)/q)			//Keep the same dimension as k in other RANS model
	)
      );	
}   //end const           	
                	

template<class BasicTurbulenceModel>
void QR<BasicTurbulenceModel>::correctNut()
{               	
    volScalarField k(this->k(fvc::grad(this->U_)));

    //this->nut_ = Ck_*sqr(this->delta())*k;		//correctNut is eddy viscosity. QR's Ck is equal to Ck^2 in Smagorinsky.C
    this->nut_ = Ck_*this->delta()*sqrt(k); 		//keep the same dimension as k in other RANS model
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/*
The template parameter of QR class is BasicTurbulenceModel. 
The template parameter of LESeddyViscosity class is BasicTurbulenceModel.
QR inherits from LESeddyViscosity base class.
The parameterized constructor of the base class can only be called using initializer list, here is an exampel
    B::B(type x):A(x){}   where the derived class is B, the base class is A.   
*/

template<class BasicTurbulenceModel>
QR<BasicTurbulenceModel>::QR			//Constructor function has no return values neither constructor prototype declaration within the class nor constructor definition.
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

    Ck_                                        //Initializing the member data/variable Ck_ of QR class in Initializer list
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
bool QR<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        Ck_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> QR<BasicTurbulenceModel>::epsilon() const
{
    volScalarField k(this->k(fvc::grad(this->U_)));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->Ce_*k*sqrt(k)/this->delta()
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> QR<BasicTurbulenceModel>::omega() const
{
    volScalarField k(this->k(fvc::grad(this->U_)));
    volScalarField epsilon(this->Ce_*k*sqrt(k)/this->delta());

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("omega", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        epsilon/(0.09*k)
    );
}


template<class BasicTurbulenceModel>
void QR<BasicTurbulenceModel>::correct()
{
    LESeddyViscosity<BasicTurbulenceModel>::correct();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
