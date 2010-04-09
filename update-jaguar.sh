#!/bin/bash
rsync -av . jaguar:/ccs/home/wirawan0/pyqmc/. -f '+ *.py' -f '+ /*/' -f '- *' 
