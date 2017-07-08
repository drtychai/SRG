"""
Author: Justin Angra
Advisor: Richard Funstahl
Last Updated: October 10, 2013

Usage: srg.py [-h] -m KINETIC_FILE -v POTENTIAL_FILE 
"""

import sys
import argparse
import numpy as np # for time-grid
import matplotlib.pyplot as plt # plotting
from time import clock  # timeing
from scipy.integrate import odeint # solving ODEs
plt.ion()

def myComm(op1, op2):
    """Returns the commutator relation for given input.
    For commuting values function returns 0.

    USAGE:
        myComm(op1, op2)

    INPUT:
        op1    - matrix of the same dimentions as op2

        op2    - matrix of the same dimentions as op1
    OUTPUT:
             - matrix representing the commutator relation between op1 and op2
    """
    return np.dot(op1, op2)-np.dot(op2, op1)

def funct_dH(T, V):
    """
    Return Hs (Hamiltonian)
    
    USAGE:
        funct_dH(T,V)
        
    INPUT:
        T   - matrix representing the eigenvalues
        
        V   - matrix representing values of the evolving potential
    """
    V = V.reshape(len(T),len(T))
    H = T + V
    eta = myComm(T, H)
    dV = myComm(eta, H)
    return dV

def f(V, mesh, T):
    """Returns dH/ds with initial values plugged in"""
    dV = funct_dH(T,V)
    arr_dV = np.array(dV).flatten()
    return arr_dV 
                                                        
def plotEvolution(soln, mesh, n_eqs):
    """Plots the time evolution of de-coupled DEs with respect
    to mesh (time-grid)
    
    USAGE:
        plotData(soln, mesh)

    INPUT:

        soln  - array of array of numerical solutions
        
        mesh  - time grid to plot over
        
        n_eqs - number of equations to plot 
    """
    
    # Plot diagonal and off diagonal elements of soln[:] 
    plt.figure()
    for i in xrange(n_eqs**2):
        plt.plot(mesh, soln[:,i], label='H%s' % i) # plot with corresponding name
    plt.axhline(y=0, color ='black', linestyle = 'dotted') # y=0 line
    plt.xlabel('lambda = 1/s ^(1/4)') # x-axis label
    plt.ylabel('Energy') # y-axis label
    plt.title('SRG Evolution') # title
    plt.legend(loc=(.6,.5), ncol=2, shadow=True) # location of legend
    plt.show()
    return

def main():
    # parser to get arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--matrix', dest='matrix_file', required=True)
    parser.add_argument('-v', '--potential', dest='vector_file', required=True)
    globals().update(vars(parser.parse_args()))

    # Read Matrix T and V
    try:
        T = np.loadtxt(matrix_file, dtype=np.float)  # Read T matrix file
        V = np.loadtxt(vector_file, dtype=np.float)  # Read V matrix file
    except Exception as e:
        print e
        sys.exit(2)
        
    # Check if both matrices are of the same size
    assert(len(T) == len(V)), "Matrices do not have same dimensions" 
    
    # Check that V is symmetric
    np.testing.assert_equal(V, np.transpose(V), err_msg = "\nPotential Matrix is not symmetric\n")

    # Calculate Matrix Size
    mat_size = len(T)
    
    v0 = np.array(V).flatten() # convert matrix to array 
    
    sMesh  = np.linspace(0, 1, 10e3)   # time grid
   
    # Solve coupled DE's    
    soln = odeint(f, v0, sMesh, args=(T,)) # returns dictionary of solutions
    
    # Create lambda mesh to take s->infinity
    lambMesh = map(lambda x: 1./x**(1./4), sMesh)
    
    # Plot functions with respect to lambda
    plotEvolution(soln, lambMesh, mat_size)
    
    print "\nSuccess"
    
if __name__ == "__main__":
    start = clock()
    main()
    print "Total run time: {0} seconds".format(clock() - start)
