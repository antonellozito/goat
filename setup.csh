#! /bin/tcsh -f

# Preliminary setup script for running goat in standalone mode (for 
# running as a submodule, see the overarching module setup workflow)
# Very rudimentary and simple, but should get the job done. 
# The following actions are done:
# - the path to 'goatrun' is added to the $PATH environment variable. 
# 'goatrun' is the main running program.
# - several environment variables are set to default values

# Say hi to the user
echo ' Welcome to the ... '
echo '    ______     _     __                                                   '
echo '   / ____/____(_)___/ /                                                   '
echo '  / / __/ ___/ / __  /                                                    '
echo ' / /_/ / /  / / /_/ /                                                     ' 
echo ' \____/_/  /_/\__,_/                                                      '
echo '    ____        __  _           _             __  _                ___    '
echo '   / __ \____  / /_(_)___ ___  (_)___  ____ _/ /_(_)___  ____     ( _ )   '
echo '  / / / / __ \/ __/ / __ `__ \/ /_  / / __ `/ __/ / __ \/ __ \   / __ \/| '
echo ' / /_/ / /_/ / /_/ / / / / / / / / /_/ /_/ / /_/ / /_/ / / / /  / /_/  <  '
echo ' \____/ .___/\__/_/_/ /_/ /_/_/ /___/\__,_/\__/_/\____/_/ /_/   \____/\/  '
echo '     /_/       __            __        __  _                              '
echo '    /   | ____/ /___ _____  / /_____ _/ /_(_)___  ____                    '
echo '   / /| |/ __  / __ `/ __ \/ __/ __ `/ __/ / __ \/ __ \                   '
echo '  / ___ / /_/ / /_/ / /_/ / /_/ /_/ / /_/ / /_/ / / / /                   '                 
echo ' /_/  |_\__,_/\__,_/ .___/\__/\__,_/\__/_/\____/_/ /_/                    '                
echo '   ______         /_/____                                                 '               
echo '  /_  __/___  ____  / / /_  ____  _  __                                   '              
echo '   / / / __ \/ __ \/ / __ \/ __ \| |/_/                                   '             
echo '  / / / /_/ / /_/ / / /_/ / /_/ />  <                                     '            
echo ' /_/  \____/\____/_/_.___/\____/_/|_|                                     '
echo ' '
echo 'Documentation: currently none!'
echo 'Current implementation: only grid optimization'

# Set goat top directory
setenv LAST_COMMAND `echo $_`
if (`echo ${LAST_COMMAND}` == "") then
  setenv GOATTOP $PWD
else
  setenv SETUP_FILE `echo ${LAST_COMMAND} | cut -d " " -f 2`
  setenv REAL_FILE `eval echo ${SETUP_FILE}`
  setenv SETUP_PATH `dirname ${REAL_FILE}`
  setenv GOATTOP `cd ${SETUP_PATH}; pwd -L`
endif

echo 'Running at goat top directory: '  ${GOATTOP}

# Add scripts to the path variable
setenv GOAT_SCRIPTPATHS ${GOATTOP}/scripts
setenv GOAT_EXECUTABLES ${GOATTOP}/executables
setenv PATH  "${GOAT_SCRIPTPATHS}:${GOAT_EXECUTABLES}:${PATH}"

# Set default environment variables
setenv COMPILER gfortran

# Set some aliases
alias gtop cd ${GOATTOP}