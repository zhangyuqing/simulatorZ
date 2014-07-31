#!/bin/sh

cd ..
R --vanilla <<RSCRIPT
library(inlinedocs);
package.skeleton.dx("zmatrix", excludePattern="AllClasses|nonexports|sn")
RSCRIPT

## Note: After doing ./runinline.sh, please remove the line \alias{doppelgangR} from man/doppelgangR-package.Rd.
grep -v '\\alias{Z.simulator}' Z.simulator/man/Z.simulator-package.Rd > Z.simulator/man/Z.simulator-package.Rd_tmp
mv Z.simulator/man/Z.simulator-package.Rd_tmp Z.simulator/man/Z.simulator-package.Rd

R CMD build Z.simulator

