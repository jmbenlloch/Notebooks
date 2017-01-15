#!/usrbin/env bash

COMMAND=$1
PYTHON_VERSION=$2

echo ${PYTHON_VERSION}


function python_version_env {
    # Activate the relevant conda env
    source activate IC${PYTHON_VERSION}
    # Set NB environment variables and download database
    nb_env
}


function nb_env {
	if [ -z $ICTDIR ]; then
		echo "ERROR: $ICTDIR is not set, run your IC manage.sh first"
	else
    	export NBDIR=`pwd`
	    echo NBDIR set to $NBDIR
    
    	echo setting PYTHONPATH
	    export PYTHONPATH=$ICTDIR
	    echo PYTHONPATH is set to $PYTHONPATH
	fi
}



## Main command dispatcher

case $COMMAND in
    python_version_env)      python_version_env ;;
    nb_install)              nb_install ;;
    
    *) echo Unrecognized command: ${COMMAND}
       echo
       echo Usage:
       echo
       echo source $0 python_version_env X.Y
       echo soruce $0 nb_env
       ;;
esac
