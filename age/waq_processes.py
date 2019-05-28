import os

from stompy.model.delft import nefis
from ctypes import * 
from ctypes.util import find_library

##

find_library('nefis')


libdir='/home/rusty/src/dfm/1.5.0/lib'

os.environ['LD_LIBRARY_PATH']=libdir

intlc=cdll.LoadLibrary(os.path.join(libdir,"libintlc.so.5")
imf=cdll.LoadLibrary("/home/rusty/src/dfm/1.5.0/lib/libimf.so")
irng=cdll.LoadLibrary("/home/rusty/src/dfm/1.5.0/lib/libirng.so")
dll=cdll.LoadLibrary('/home/rusty/src/dfm/1.5.0/lib/libnefis.so.0.0.0')
