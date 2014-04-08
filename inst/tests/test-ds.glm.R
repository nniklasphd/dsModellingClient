#-------------------------------------------------------------------------------
# Copyright (c) 2014 OBiBa. All rights reserved.
#  
# This program and the accompanying materials
# are made available under the terms of the GNU Public License v3.0.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------

#
# Set up
#

context("dsbaseclient::ds.glm")

options(datashield.variables=list("DIS_DIAB","PM_BMI_CONTINUOUS","LAB_HDL", "GENDER"))
source("setup.R")

#
# Tests
#

context("dsbaseclient::ds.glm() run a GLM without interaction (e.g. diabetes prediction using BMI and HDL levels and GENDER)")
mod <- ds.glm(datasources=opals,formula=D$DIS_DIAB~D$PM_BMI_CONTINUOUS+D$LAB_HDL+D$GENDER,family=quote(binomial))
# TODO do more than a smoke test

context("dsbaseclient::ds.glm()  run the above GLM model with an intercept (eg. intercept = 1)")
mod <- ds.glm(datasources=opals,formula=D$DIS_DIAB~1+D$PM_BMI_CONTINUOUS+D$LAB_HDL+D$GENDER,family=quote(binomial))
# TODO do more than a smoke test

context("dsbaseclient::ds.glm() run the above GLM with interaction HDL and GENDER")
mod <- ds.glm(datasources=opals,formula=D$DIS_DIAB~D$PM_BMI_CONTINUOUS+D$LAB_HDL*D$GENDER,family=quote(binomial))
# TODO do more than a smoke test

context("dsbaseclient::ds.glm() now run the same GLM but with interaction between BMI and HDL")
mod <- ds.glm(datasources=opals,formula=D$DIS_DIAB~D$PM_BMI_CONTINUOUS*D$LAB_HDL+D$GENDER,family=quote(binomial))
# TODO do more than a smoke test

#
# Tear down
#

source("teardown.R")