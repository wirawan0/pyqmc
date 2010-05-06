#!/bin/bash
rsync -av . jaguar-ext:/ccs/home/wirawan0/pyqmc/. -f '+ *.py' -f '+ *.basis' -f '+ /**/' -f '- *'
