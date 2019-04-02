#-------------------------------------------------------------------------------
# Copyright (c) 2019 OBiBa. All rights reserved.
#  
# This program and the accompanying materials
# are made available under the terms of the GNU Public License v3.0.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------

#
# Datashield test suite set up
#

library(dsModellingClient)
library(DSLite)
library(testthat)

# load test datasets, instanciate a DSLiteServer and corresponding logindata in current environment
logindata <- DSLite::setupCNSIMTest("dsBase", env = environment())

myvar <- list("LAB_TSC", "LAB_HDL")
conns <- datashield.login(logins=logindata,assign=TRUE,variables=getOption("datashield.variables", myvar))

