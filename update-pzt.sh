#!/bin/bash
rsync -av . pmn:/home/wirawan0/local/lib/python/pyqmc/. -f '+ *.py' -f '+ *.basis' -f '+ /**/' -f '- *'
