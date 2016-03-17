###################################################
# START BEAGLE environment variables
###################################################
setenv BEAGLE_ROOT "${HOME}/beagle/files"

setenv BEAGLE_PARAM_DIR "${BEAGLE_ROOT}/param"
setenv BEAGLE_TEST_FILES "${BEAGLE_ROOT}/tests"

setenv BEAGLE_TEMPLATES "${BEAGLE_ROOT}/templates"
setenv BEAGLE_FILTERS "${BEAGLE_ROOT}/filters"
setenv BEAGLE_DUST "${BEAGLE_ROOT}/dust"
setenv BEAGLE_DATA "${BEAGLE_ROOT}/data"
setenv BEAGLE_SF_CHE "${BEAGLE_ROOT}/sf_che"
setenv BEAGLE_RESULTS "${BEAGLE_ROOT}/results"

# This environment variable is already defined by BC03 routines, but here you
# redefine it to be sure that the filters used are those in the BEAGLE directory
# tree 
setenv FILTERS "${BEAGLE_FILTERS}/FILTERBIN.RES"

###################################################
# END BEAGLE environment variables
###################################################
