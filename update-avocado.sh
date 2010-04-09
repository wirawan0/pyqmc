#!/bin/bash
rsync -av . avo:/home/wirawan0/lib/python/pyqmc/. -f '+ *.py' -f '+ *.basis' -f '+ /**/' -f '- *' 
