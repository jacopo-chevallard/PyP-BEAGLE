###################################################
# START BEAGLE environment variables
###################################################
export BEAGLE_ROOT="${HOME}/beagle/files"

export BEAGLE_PARAM_DIR="${BEAGLE_ROOT}/param"
export BEAGLE_TEST_FILES="${BEAGLE_ROOT}/tests"

export BEAGLE_TEMPLATES="${BEAGLE_ROOT}/templates"
export BEAGLE_FILTERS="${BEAGLE_ROOT}/filters"
export BEAGLE_DUST="${BEAGLE_ROOT}/dust"
export BEAGLE_DATA="${BEAGLE_ROOT}/data"
export BEAGLE_SF_CHE="${BEAGLE_ROOT}/sf_che"
export BEAGLE_RESULTS="${BEAGLE_ROOT}/results"

# This environment variable is already defined by BC03 routines, but here you
# redefine it to be sure that the filters used are those in the BEAGLE directory
# tree 
export FILTERS="${BEAGLE_FILTERS}/FILTERBIN.RES"

###################################################
# END BEAGLE environment variables
###################################################
