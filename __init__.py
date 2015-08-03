sub_paths =  ["", "sedov", "analysis"]

import os
__path__ = [os.path.join(__path__[0], path) for path in sub_paths]
__path__ = [os.path.normpath(path) for path in __path__]
