from . import fhelmholtz
import numpy as np

# list of common blocks that hold interesting information
CBLOCK_NAMES = ('crpc1', 'deedoo', 'etotc1', 'etapc1', 'ptotc1',
                'stotc1', 'th_xni_ion', 'thcou', 'thdergc1', 'thdergc2',
                'thdertc1', 'thdertc2', 'theepc1', 'thegasc2', 'theion',
                'therad', 'thetaion', 'thinp', 'thmax', 'thpepc1',
                'thpgasc1', 'thpion', 'thprad', 'thsepc1', 'thsgasc1',
                'thsion', 'thsrad', 'thxip', 'thxnec1')


class HelmholtzOutput:

    def __init__(self, size, shape):

        # set size and shape data
        self.size = size
        self.shape = shape

        # loop through and nicely reformat everything
        for cblock_name in CBLOCK_NAMES:
            cblock = getattr(fhelmholtz,cblock_name)
            for row_name in vars(cblock):
                row_data = np.copy(getattr(cblock,row_name))
                setattr(self, self._demangle(row_name), self._reshape(row_data))

    def _demangle(self, name):
        # remove the "_row" postfix
        return name[:-4]

    def _reshape(self, data):
        # put things back like they came in
        return np.reshape(data[:self.size], self.shape)


def _make_uniform_arrays(inputs):

    # make numpy arrays out of everthing
    arrays = [np.array(array) for array in inputs]

    # set size & shape to 1st non-scalar input
    size = 1
    shape = (1,)
    for array in arrays:
        if array.size != 1:
            size  = array.size
            shape = array.shape
            break

    outputs = []
    for array in arrays:
        if array.size  == 1:
            outputs.append(np.tile(array.flatten(), size))
        else:
            outputs.append(array.flatten())

    return size, shape, outputs


def helmeos(dens, temp, abar, zbar):

    # make sure everything is the same size and shape
    inputs = (dens, temp, abar, zbar)
    size, shape, finputs = _make_uniform_arrays(inputs)

    # call the eos
    fhelmholtz.call_helmeos(*finputs)

    # container for output
    return HelmholtzOutput(size, shape)


def helmeos_DE(dens, ener, abar, zbar, tguess = None):

    # set default temperature guess
    if tguess == None:
        tguess = 1e7

    # make sure everything is the same size and shape
    inputs = (dens, ener, abar, zbar, tguess)
    size, shape, finputs = _make_uniform_arrays(inputs)

    # call the eos
    fhelmholtz.call_helmeos_de(*finputs)

    return HelmholtzOutput(size, shape)


def helmeos_DP(dens, pres, abar, zbar, tguess = None):

    # set default temperature guess
    if tguess is None:
        tguess = 1e7

    # make sure everything is the same size and shape
    inputs = (dens, pres, abar, zbar, tguess)
    size, shape, finputs = _make_uniform_arrays(inputs)

    # call the eos
    fhelmholtz.call_helmeos_dp(*finputs)

    return HelmholtzOutput(size, shape)


def helmeos_DS(dens, entr, abar, zbar, tguess = None):

    # set default temperature guess
    if tguess is None:
        tguess = 1e7

    # make sure everything is the same size and shape
    inputs = (dens, entr, abar, zbar, tguess)
    size, shape, finputs = _make_uniform_arrays(inputs)

    # call the eos
    fhelmholtz.call_helmeos_ds(*finputs)

    return HelmholtzOutput(size, shape)
