#!/bin/sh

cd ..
R --vanilla <<RSCRIPT
library(inlinedocs);
package.skeleton.dx("simulatorZ", excludePattern="AllClasses|nonexports|sn")
RSCRIPT

## Note: After doing ./runinline.sh, please remove the line \alias{doppelgangR} from man/doppelgangR-package.Rd.
grep -v '\\alias{simulatorZ}' simulatorZ/man/simulatorZ-package.Rd > simulatorZ/man/simulatorZ-package.Rd_tmp
mv simulatorZ/man/simulatorZ-package.Rd_tmp simulatorZ/man/simulatorZ-package.Rd

R CMD build simulatorZ

