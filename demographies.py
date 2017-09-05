"""
Custom demographic models.
"""
import numpy
import dadi
import sys

def trunk_2epoch_sizechange((nu_anc_s, TB), n1, pts):
    """
    Model with single size change in single population 

    ~~parameters to maximize~~
    nu_anc_s: size change strength
    TB: time between size change and present

    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    
    #size change to size nu_anc_S, remains this size for time period TB
    phi=dadi.Integration.one_pop(phi,xx,TB,nu=nu_anc_s)
    
    #calculate the spectrum.
    sfs = dadi.Spectrum.from_phi(phi, n1, (xx,))
    
    return sfs

def trunk_3epoch_sizechange((nu_anc_s, TB,TR), n1, pts):
    """
    Model with 2 size changes in single population
    
    ~~parameters to maximize~~    
    nu_anc_s: size change strength
    TB: duration of size change
    TR: time between recovery and present

    ~~fixed parameters~~
    nuR: size of population after recovery
    
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    
    #size change to size nuS, remains for time period TB
    phi=dadi.Integration.one_pop(phi,xx,TB,nu=nu_anc_s)
    nuR=1
    phi=dadi.Integration.one_pop(phi,xx,TR,nu=nuR)
    
    #calculate the spectrum.
    sfs = dadi.Spectrum.from_phi(phi, n1, (xx,))
    
    return sfs

def bneck_3params((nu2B, TB, TF), (n1,n2), pts):
    """
    Model with split, bottleneck in pop2, recovers to 1
        
    ~~parameters to maximize~~    
    nu2B: The bottleneck size for pop2
    TB: duration of bottleneck
    TF: The time between the bottleneck recovery and present

    ~~fixed parameters~~
    nu1: population size of population 1, also the size of population 2 after bottleneck recovery
        
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx =yy= dadi.Numerics.default_grid(pts)
    
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    #immediate bottleneck
    nu1=1
    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nu1, nu2=nu2B,m12=0,m21=0)
    phi = dadi.Integration.two_pops(phi, xx, TF, nu1=nu1, nu2=nu1,m12=0,m21=0)
    
    # Finally, calculate the spectrum.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return sfs

def bneck_2params((nu2B, TF), (n1,n2), pts):
    """
    Model with split, bottleneck in pop2, recovers to 1, length of bottleneck is fixed
    
    ~~parameters to maximize~~ 
    nu2B: The bottleneck strength for pop2
    TF: The time between the recovery and present

    ~~fixed parameters~~
    TB: length of bottleneck
    nu1: population size of population 1, also the size of population 2 after bottleneck recovery
    
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx =yy= dadi.Numerics.default_grid(pts)
    
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    #immediate bottleneck
    TB=4200/(2.0*50000)
    nu1=1
    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nu1, nu2=nu2B,m12=0,m21=0)
    phi = dadi.Integration.two_pops(phi, xx, TF, nu1=nu1, nu2=nu1,m12=0,m21=0)
    
    # Finally, calculate the spectrum.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return sfs

def trunk_2epoch_sizechange_bneck_3param((nu_anc_s, t_anc_s, TB,nu2B,TF), (n1,n2), pts):
    """
    Model with single size change in ancestral population, divergence, 3 param bottleneck in pop 2
    
    ~~parameters to maximize~~    
    nu_anc_s: size change in ancestral population
    t_anc_s: time between ancestral size change and population split
    TB: duration of bottleneck in population 2
    nu2B: strength of bottleneck in pop2
    TF: time between pop 2 recovery and present

    ~~fixed parameters~~
    nu1: size of pop1 1 and size of pop2 after recovery from bottleneck. Fixed at size nu_anc_s which is size of ancestral population after size change
    
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    # Define the grid we'll use
    xx=yy = dadi.Numerics.default_grid(pts)
    
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    
    #size change to size nu_anc_s, remains for time period t_anc_s
    phi=dadi.Integration.one_pop(phi,xx,t_anc_s,nu=nu_anc_s)
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    #immediate bottleneck
    nu1=nu_anc_s 
    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nu1, nu2=nu2B,m12=0,m21=0)
    phi = dadi.Integration.two_pops(phi, xx, TF, nu1=nu1, nu2=nu1,m12=0,m21=0)
    
    # Finally, calculate the spectrum.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return sfs

def trunk_2epoch_sizechange_bneck_2param((nu_anc_s, t_anc_s,nu2B,TF), (n1,n2), pts):
    """
    Model with single size change in ancestral population, divergence, bottleneck in pop 2, length of bottleneck is fixed
    
    ~~parameters to maximize~~
    nu_anc_s: size change in ancestral population
    t_anc_s: time between ancestral size change and population split
    nu2B: strength of bottleneck in pop2
    TF: time between pop 2 recovery and present

    ~~fixed parameters~~
    nu1: size of pop1 1 and size of pop2 after recovery from bottleneck. Fixed at size nu_anc_s which is size of ancestral population after size change
    TB: length of bottleneck in population 2
    
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx =yy= dadi.Numerics.default_grid(pts)
    
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    
    #size change to size nu_anc_s, remains for time period t_anc_s
    phi=dadi.Integration.one_pop(phi,xx,t_anc_s,nu=nu_anc_s)
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    #immediate bottleneck
    TB=4200/(2.0*50000)
    nu1=nu_anc_s
    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nu1, nu2=nu2B,m12=0,m21=0)
    phi = dadi.Integration.two_pops(phi, xx, TF, nu1=nu1, nu2=nu1,m12=0,m21=0)
    
    # Finally, calculate the spectrum.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return sfs

def trunk_3epoch_sizechange_bneck_3param((nu_anc_s, t_anc_s,t_anc_r, TB,nu2B,TF), (n1,n2), pts):
    """
    Model with single size change in ancestral population, divergence, bottleneck in pop 2
    
    ~~parameters to maximize~~
    nu_anc_s: size change in ancestral population
    t_anc_s: duration of size change in ancestral population
    t_anc_r: time between recovery and population split 
    TB: duration of bottleneck in population 2
    nu2B: strength of bottleneck in pop2
    TF: time between pop 2 recovery and present

    ~~fixed parameters~~
    nu1: size of pop1 1 and size of pop2 after recovery from bottleneck
    
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx =yy= dadi.Numerics.default_grid(pts)
    
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    
    #size change to size nu_anc_s, remains for time period t_anc_s
    phi=dadi.Integration.one_pop(phi,xx,t_anc_s,nu=nu_anc_s)
    nu1=1
    phi=dadi.Integration.one_pop(phi,xx,t_anc_r,nu=nu1)
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    #immediate bottleneck
    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nu1, nu2=nu2B,m12=0,m21=0)
    phi = dadi.Integration.two_pops(phi, xx, TF, nu1=nu1, nu2=nu1,m12=0,m21=0)
    
    # Finally, calculate the spectrum.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return sfs

def trunk_3epoch_sizechange_bneck_2param((nu_anc_s, t_anc_s,t_anc_r,nu2B,TF), (n1,n2), pts):
    """
    Model with single size change in ancestral population, divergence, bottleneck in pop 2
    
    ~~parameters to maximize~~
    nu_anc_s: size change in ancestral population
    t_anc_s: duration of size change in ancestral population
    t_anc_r: time between recovery and population split
    nu2B: strength of bottleneck in pop2
    TF: time between pop 2 recovery and present

    ~~fixed parameters~~
    nu1: size of pop1 1 and size of pop2 after recovery from bottleneck
    TB: length of bottleneck in population 2
    
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx =yy= dadi.Numerics.default_grid(pts)
    
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    
    #size change to size nu_anc_s, remains for time period t_anc_s
    phi=dadi.Integration.one_pop(phi,xx,t_anc_s,nu=nu_anc_s)
    nu1=1
    phi=dadi.Integration.one_pop(phi,xx,t_anc_r,nu=nu1)
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    #immediate bottleneck
    TB=4200/(2.0*50000)
    phi = dadi.Integration.two_pops(phi, xx, TB, nu1=nu1, nu2=nu2B,m12=0,m21=0)
    phi = dadi.Integration.two_pops(phi, xx, TF, nu1=nu1, nu2=nu1,m12=0,m21=0)
    
    # Finally, calculate the spectrum.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return sfs

def IM_2params((T,m),(n1,n2),pts):
    """
    Model with split into equal population sizes with no size change, symmetric migration

    ~~parameters to maximize~~
    T: time since split
    m: migration rate

    ~~fixed parameters~~
    nu: population sizes after split, no size change so it remains 1

    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """

    # Define the grid we'll use
    xx =yy= dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    nu=1
    phi = dadi.Integration.two_pops(phi,xx,T,nu1=nu,nu2=nu,m12=m,m21=m)
    
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return sfs


def trunk_2epoch_sizechange_IM_2param((nu_anc_s, t_anc_s,T,m), (n1,n2), pts):
    """
    Model with population size change in ancestral population, split into equal population sizes with no size change, symmetric migration

    ~~parameters to maximize~~
    nu_anc_s: size change in ancestral population
    t_anc_s: time between ancestral size change and population split
    T: time since split
    m: migration rate

    ~~fixed parameters~~
    nu: population sizes after split

    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx =yy= dadi.Numerics.default_grid(pts)
    
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    
    #size change to size nu_anc_s, remains for time period t_anc_s
    phi=dadi.Integration.one_pop(phi,xx,t_anc_s,nu=nu_anc_s)
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    nu=nu_anc_s
    phi = dadi.Integration.two_pops(phi,xx,T,nu1=nu,nu2=nu,m12=m,m21=m)
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return sfs

def trunk_3epoch_sizechange_IM_2param((nu_anc_s, t_anc_s,t_anc_r,T,m), (n1,n2), pts):
    """
    Model with population size change and recovery in ancestral population, split into equal population sizes with no size change, symmetric migration

    ~~parameters to maximize~~
    nu_anc_s: size change in ancestral population
    t_anc_s: duration of size change in ancestral population
    t_anc_r: time between recovery and population split
    T: time since split
    m: migration rate

    ~~fixed parameters~~
    nu: population sizes after split

    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx =yy= dadi.Numerics.default_grid(pts)
    
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    
    #size change to size nu_anc_s, remains for time period t_anc_s
    phi=dadi.Integration.one_pop(phi,xx,t_anc_s,nu=nu_anc_s)
    nu1=1
    phi=dadi.Integration.one_pop(phi,xx,t_anc_r,nu=nu1)
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    nu=1
    phi = dadi.Integration.two_pops(phi,xx,T,nu1=nu,nu2=nu,m12=m,m21=m)
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return sfs