
#define PY_SSIZE_T_CLEAN
#include <Python.h>
/*  Based on the description of the EEP 3.x file format in:
 *  * cnt_riff.txt by Rainer Nowagk & Maren Grigutsch.
 */


/* Start of old raweep.h header file */

#include <stdint.h>
#include <stdio.h>

#if defined ( __GNUC__ )
    uint64_t swap_uint64 ( uint64_t val ) {
        return __builtin_bswap64 ( val );
    }
#elif defined ( _MSC_VER )
    uint64_t swap_uint64 ( uint64_t val ) {
        return _byteswap_uint64 ( val );
    }
#else
    /* By chmike in https://stackoverflow.com/a/2637138 */
    uint64_t swap_uint64 ( uint64_t val ) {
        val = ((val << 8) & 0xFF00FF00FF00FF00ULL ) | ((val >> 8) & 0x00FF00FF00FF00FFULL );
        val = ((val << 16) & 0xFFFF0000FFFF0000ULL ) | ((val >> 16) & 0x0000FFFF0000FFFFULL );
        return (val << 32) | (val >> 32);
    }
#endif

void display ( uint32_t num ) {
    int i;
    printf ( "Value: %i\n", num );
    printf ( "Binary expression: ");
    for ( i = 0; i < 32; i++ ) {
        printf ( "%i", ( num & ( 1 << (31 - i) ) ) != 0 );
    }
    printf ( "\n" );
    printf ( "\n" );
}

void display64 ( uint64_t num ) {
    int i;
    printf ( "Value: %zi\n", num );
    printf ( "Binary expression: ");
    for ( i = 0; i < 64; i++ ) {
        printf ( "%i", ( num & ( (uint64_t)1 << (63 - i) ) ) != 0 );
    }
    printf ( "\n" );
    printf ( "\n" );
}



int32_t bin2int32 ( const uint8_t *byte, int64_t offset, uint32_t length ) {
    
    uint64_t tmpval;
    int32_t intval;
    uint32_t mask;
    
    
    /* Calculates the mask. */
    mask   = (uint32_t) ( ( (uint64_t) 1 ) << length ) - 1;
    
    /* Gets the continuous stream of bits. */
    memcpy ( &tmpval, byte + offset / 8, sizeof ( tmpval ) );
    
    /* Flips the bytes. */
    tmpval = swap_uint64 ( tmpval );
    
    /* Selects the desired chunk. */
    intval = ( tmpval >> ( 8 * sizeof ( tmpval ) - length - offset % 8 ) ) & mask;
    
    /* If the value is negative fills the leftmost zeros. */
    if ( ( intval >> ( length - 1 ) ) != 0 ) {
        intval = intval | ~mask;
    }
    
    return intval;
}


int32_t bin2uint32 ( const uint8_t *byte, int64_t offset, uint32_t length ) {
    
    uint64_t tmpval;
    uint32_t intval;
    uint32_t mask;
    
    
    /* Calculates the mask. */
    mask   = (uint32_t) ( ( (uint64_t) 1 ) << length ) - 1;
    
    /* Gets the continuous stream of bits. */
    memcpy ( &tmpval, byte + offset / 8, sizeof ( tmpval ) );
    
    /* Flips the bytes. */
    tmpval = swap_uint64 ( tmpval );
    
    /* Selects the desired chunk. */
    intval = ( tmpval >> ( 8 * sizeof ( tmpval ) - length - offset % 8 ) ) & mask;
    
    return intval;
}


void apply_residuals ( int32_t *data, uint64_t nsamp, uint32_t method ) {
    
    uint64_t index;
    
    
    /* Extendes the residuals depending onf the compression method. */
    switch ( method & 7 ) {
        
        /* For method 0 (or 8) does nothing. */
        case 0:
            break;
            
        /* For method 1 (or 9) applies the residuals along samples. */
        case 1:
            
            /* Goes through each sample and accumulates one residual. */
            for ( index = 1; index < nsamp; index++ ) {
                data [ index ] = data [ index ] + data [ index - 1 ];
            }
            
            break;
            
        /* For method 2 (or 10) applies the residuals along samples twice. */
        case 2:
            
            /* The second sample only has one residual. */
            data [1] = data [0] + data [1];
            
            /* For the rest of samples accumulates two residuals. */
            for ( index = 2; index < nsamp; index++ ) {
                data [ index ] = data [ index ] + data [ index - 1 ] + ( data [ index - 1 ] - data [ index - 2 ] );
            }
            
            break;
            
        /* For method 3 (or 11) applies the residuals along samples and channels. */
        case 3:
            
            /* Goes through each sample from the second on. */
            for ( index = 1; index < nsamp; index ++ ) {
                
                /* Accumulates one residual. */
                data [ index ] = data [ index ] + data [ index - 1 ];
                
                /* Also adds the residual of the previous channel. */
                data [ index ] = data [ index ] + ( data [ index - nsamp ] - data [ index - nsamp - 1 ] );
            }
            
            break;
    }
}


int64_t read_channel ( int32_t *data, const uint8_t *byte, int64_t off, uint64_t nsamp, uint32_t method ) {
    
    uint32_t mbit, dbit, nbit, xbit;
    uint64_t index;
    int32_t datum;
    
    
    /* Gets the length of the data and metadata depending on the method. */
    mbit   =  4 + ( ( method & 8 ) >> 2 );
    dbit   = 16 + ( ( method & 8 ) << 1 );
    
    /* Checks if the data is stored as residuals. */
    if ( ( method & 7 ) == 0 ) {
        
        /* Data is stored as samples. */
        nbit   = dbit;
        xbit   = 0;
        off   += 4;
        
    } else {
        
        /* Reads the length of the standard and extended residuals. */
        nbit   = bin2uint32 ( byte, off, mbit );
        off   += mbit;
        xbit   = bin2uint32 ( byte, off, mbit );
        off   += mbit;
        
        /* If xbit is equal to nbit there are no extended residuals. */
        if ( xbit == nbit ) xbit = 0;
    }
    
    
    /* Reads the first sample of the channel. */
    datum  = bin2int32 ( byte, off, dbit );
    off   += dbit;
    
    /* Stores the first sample for this channel. */
    data [0] = datum;
    
    
    /* Goes through the samples. */
    for ( index = 1; index < nsamp; index ++ ) {
        
        /* Reads nbit bits. */
        datum  = bin2int32 ( byte, off, nbit );
        off   += nbit;
        
        /* Checks if we need the extended residual. */
        if ( xbit && ( datum == -pow ( 2, nbit - 1 ) ) ) {
            
            /* Reads nbit bits. */
            datum  = bin2int32 ( byte, off, xbit );
            off   += xbit;
        }
        
        /* Stores the current residual. */
        data [ index ] = datum;
    }
    
    /* Aligns the offset with the byte. */
    off    = 8 * ( ( off - 1 ) / 8 + 1 );
    
    return off;
}


int64_t read_block ( int32_t *data, const uint8_t *byte, int64_t off, uint64_t nsamp, uint64_t nchan ) {
    
    uint64_t index;
    uint32_t method;
    
    
    /* Iterates through channels. */
    for ( index = 0; index < nchan; index ++ ) {
        
        /* Checks the compression method. */
        method = bin2uint32 ( byte, off, 4 );
        off   += 4;
        
        /* Only methods 0, 1, 2, 3, 8, 9, 10, 11 are allowed. */
        if ( ( method & 7 ) > 3 )
            return -1;
        
        /* The first channel cannot have method 3/11. */
        if ( index == 0 && ( ( method & 7 ) == 3 ) )
            return -2;
        
        /* Reads the data for this channel. */
        off    = read_channel ( data + ( index * nsamp ), byte, off, nsamp, method );
        
        /* Expands the residuals for this channel. */
        apply_residuals ( data + ( index * nsamp ), nsamp, method );
    }
    
    return off;
}

/* End of old raweep.h header file */


/*
#include "raweep.h"
*/


static PyObject *
raweep_read_block ( PyObject *self, PyObject *args ) {

    const uint8_t * bytes;
    const size_t bsize;
    const uint64_t nsamp, nchan, offset;
    int64_t off;
    int32_t * data;
    PyObject * result;

    if ( !PyArg_ParseTuple ( args, "y#KKK", &bytes, &bsize, &nsamp, &nchan, &offset ) )
        return NULL;
    
    
    /* Reserves memory for the block data. */
    data = PyMem_Malloc ( nsamp * nchan * sizeof ( *data ) );
    
    /* Reads the block. */
    off  = read_block ( data, bytes, offset, nsamp, nchan );
    
    if ( off < 1 )
        PyErr_SetString ( PyExc_RuntimeError, "The compression method for the first channel cannot be inter-channel residuals (method 3/11)." );
    
    
    /* Destroys the archived data. */
    /*PyMem_Free ( data );*/

    /* Returns the raw data stream and the offset. */
    /* return Py_BuildValue ( "y#K", data, nchan * nsamp * sizeof ( uint32_t ), off ); */
    /* Avoid memory leak - Victor */

    result = Py_BuildValue ( "y#K", data, nchan * nsamp * sizeof ( uint32_t ), off );
    PyMem_Free ( data );
    return result;
}



/* Defines the method table (available methods). */
static PyMethodDef
raweep_methods [] = {
    { "read_block",  raweep_read_block, METH_VARARGS, "Decodes a block of data from a raw3 stream." },
    { NULL, NULL, 0, NULL }        /* Sentinel */
};

/* Describes the module */
static struct
PyModuleDef raweepmodule = {
    PyModuleDef_HEAD_INIT,
    "raweep", 
    PyDoc_STR ( "Module to read raw data from EEP (CNT) files." ),
    -1,
    raweep_methods
};

/* Module's initialization function. */
PyMODINIT_FUNC
PyInit_raweep ( void ) {
    return PyModule_Create ( &raweepmodule );
}


int
main ( int argc, char *argv [] ) {
    wchar_t *program = Py_DecodeLocale ( argv [0], NULL );
    if ( program == NULL ) {
        fprintf ( stderr, "Fatal error: cannot decode argv[0]\n" );
        exit (1);
    }

    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName ( program );

    /* Add a built-in module, before Py_Initialize */
    if ( PyImport_AppendInittab ( "raweep", PyInit_raweep ) == -1 ) {
        fprintf ( stderr, "Error: could not extend in-built modules table\n" );
        exit (1);
    }

    /* Initializes the Python interpreter. */
    Py_Initialize ();

    /* Imports the module. */
    PyObject *pmodule = PyImport_ImportModule ( "raweep" );
    
    if ( !pmodule ) {
        PyErr_Print ();
        fprintf ( stderr, "Error: could not import module 'raweep'\n" );
    }


    PyMem_RawFree(program);
    return 0;
}
