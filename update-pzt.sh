#!/bin/bash
rsync -av . pzt:/home/wirawan0/local/lib/python/pyqmc/. -f '+ *.py' -f '+ /*/' -f '- *' 
