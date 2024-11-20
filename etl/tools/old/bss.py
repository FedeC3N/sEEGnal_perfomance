#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 18:53:36 2023

@author: bru
"""
'''
try:
    profile
except:
    import line_profiler
    profile = line_profiler.LineProfiler()
'''

import scipy
import numpy
import time

'''
@profile'''
def sobi ( data, nlag = None, nsource = None ):
    
    '''
    Based on function:
        * sobi.m by A. Belouchrani and A. Cichocki.
    
    Based on publications:
        * Belouchrani, Abed-Meraim, Cardoso & Moulines 1993
          Proc. Int. Conf. on Digital Sig. Proc. 346-351.
        * Belouchrani, Abed-Meraim 1993
          Proc. Gretsi, (Juan-les-pins) 309-312.
        * Cichocki & Amari 2003
          Adaptive Blind Signal and Image Processing, Wiley.
    '''
    
    
    # Gets the metadata (assumes data in MNE-Python shape).
    nchannel = data.shape [-2]
    nsample  = data.shape [-1]
    data     = data.reshape ( ( -1, nchannel, nsample ) )
    nepoch   = data.shape [0]
    
    
    # If no requested otherwise, sets the number of lags to 1/3 of the data.
    if nlag is None:
        nlag     = min ( 100, numpy.ceil ( nsample / 3 ) )
    
    # If no requested otherwise, uses as many sources as channels.
    if nsource is None:
        nsource  = nchannel
    
    
    # Rewrites the data as channels by samples by epochs.
    data     = numpy.moveaxis ( data, 0, -1 ).copy ()
    
    # Demeans the data.
    data     = data - data.mean ( axis = 1, keepdims = True )
    
    # Concatenates all the epochs.
    data     = data.reshape ( ( data.shape [0], -1 ) )
    
    
    # Calculates a whitener using SVD.
    [ u, s, v ] = scipy.linalg.svd ( numpy.dot ( data, data.T ) )
    
    # Restricts the number of sources, if requested.
    u        = u [ :, : nsource ]
    s        = s [ : nsource ]
    
    # Calculates the whitener and unwhitener,
    wfilter  = ( u / ( numpy.sqrt ( s ) ) ).T
    uwfilter = scipy.linalg.pinv ( wfilter )
    
    
    # Whitens the data.
    wdata    = numpy.dot ( wfilter, data )
    
    # Restores the original shape of the data.
    wdata    = wdata.reshape ( ( nchannel, nsample, nepoch ) )
    
    # Rewrites the whitened data as samples by epochs by sources.
    wdata    = numpy.moveaxis ( wdata, 0, -1 ).copy ()
    
    
    # Initializes the correlations matrix.
    xcorr    = numpy.zeros ( ( nchannel, nchannel, nlag ) )
    
    # Gets the correlation matrix for each lag.
    for lag in range ( nlag ):
        
        # Gets the overlapping data for the current lag.
        wdata1   = wdata [ lag + 1:, :, : ].reshape ( ( -1, nsource ) )
        wdata2   = wdata [ : nsample - lag - 1, :, : ].reshape ( ( -1, nsource ) )
        
        # Gets the unbiased normalized correlation.
        Rxp      = numpy.dot ( wdata1.T, wdata2 )
        Rxp      = Rxp / nepoch / ( nsample - lag - 1 )
        Rxp      = Rxp * scipy.linalg.norm ( Rxp )
        
        # Stores the normalized correlation for this lag.
        xcorr [ :, :, lag ] = Rxp
    
    
    # Initializes the flags and counters.
    iter_    = 0
    loop     = True
    
    # Defines the stop condition for the iterations.
    eps      = 1 / numpy.sqrt ( nsample ) / 100
    eps2     = eps ** 2
    
    
    # Creates the rotation mask.
    mask     = numpy.array ( [
        [ 1, 1, 0 ],
        [ 1, 1, 0 ],
        [ 0, 0, 1 ] ] )
    
    
    # Initializes the mixing matrix to the identity.
    wmixing  = numpy.eye ( nchannel )
    
    # Performs the joint diagonalization.
    start_SOBI = time.time()
    while loop:

        # Updates the iteration index and reset the flags.
        iter_    += 1
        loop     = False
        
        #print ( iter_ )
        
    
        # Goes through each pair of sources.
        for p in range ( nchannel - 1 ):
            for q in range ( p + 1, nchannel ):

                # Check the elapsed time
                current_time_SOBI = time.time()
                elapsed_time_current_iteration = (current_time_SOBI - start_SOBI)

                # Raise an error if the SOBI process takes more than 6 hours
                if elapsed_time_current_iteration > 6 * 60 * 60:
                    raise Exception('SOBI process aborted. It took more than 6 hours')
                
                # Computes the dependence between the sources.
                g1      = xcorr [ p, p, : ] - xcorr [ q, q, : ]
                g2      = xcorr [ p, q, : ] + xcorr [ q, p, : ]
                g3      = xcorr [ q, p, : ] - xcorr [ p, q, : ]
                
                g       = numpy.array ( [ g1, g2, g3 ] )
                gg      = numpy.dot ( g, g.T )
                gg      = gg * mask
                
                
                # Gets the Givens rotation from the main eigenvector.
                angles  = scipy.linalg.svd ( gg ) [0]
                angles  = angles [ :, 0 ]
                angles  = numpy.sign ( angles [0] ) * angles
                
                c       = numpy.sqrt ( 0.5 * ( 1 + angles [0] ) )
                #sr      = 0.5 * ( angles [1] - 1j * angles [2] ) / c
                #sc      = 0.5 * ( angles [1] + 1j * angles [2] ) / c
                sr      = 0.5 * ( angles [1] ) / c
                sc      = 0.5 * ( angles [1] ) / c
                
                # Performs the rotation, if required.
                if sc * sr > eps2:
                    
                    # At least one element was modified, so iterates again.
                    loop     = True
                    
                    # Rotates the correlation matrix along the columns.
                    dummyp   = xcorr [ :, p, : ].copy ()
                    dummyq   = xcorr [ :, q, : ].copy ()
                    xcorr [ :, p, : ] = c * dummyp + sr * dummyq
                    xcorr [ :, q, : ] = c * dummyq - sc * dummyp
                    
                    # Rotates the correlation matrix along the rows.
                    dummyp   = xcorr [ p, :, : ].copy ()
                    dummyq   = xcorr [ q, :, : ].copy ()
                    xcorr [ p, :, : ] = c * dummyp + sc * dummyq
                    xcorr [ q, :, : ] = c * dummyq - sr * dummyp
                    
                    # Rotates the mixing matrix.
                    dummyp   = wmixing [ p, : ].copy ()
                    dummyq   = wmixing [ q, : ].copy ()
                    wmixing [ p, : ] = c * dummyp + sr * dummyq
                    wmixing [ q, : ] = c * dummyq - sc * dummyp





    # Un-whitens the mixing matrix.
    mixing   = numpy.dot ( uwfilter, wmixing.T );
    
    # Estimates the unmixing matrix.
    unmixing = scipy.linalg.pinv ( mixing )
    
    
    # Returns the mixing and unmixing matrices.
    return mixing, unmixing
