#!/usrbin/env bash

COMMAND=$1
PYTHON_VERSION=$2

echo ${PYTHON_VERSION}


function setup_git {
    git config --add include.path `pwd`/.gitconfig
}


function python_version_env {
    # Activate the relevant conda env
    source activate IC${PYTHON_VERSION}
    # Set NB environment variables and download database
    notebook_env
}


function notebook_env {
	if [ -z $ICTDIR ]; then
		echo "ERROR: $ICTDIR is not set, run your IC manage.sh first"
	else
        export NBDIR=`pwd`
    	export IC_NOTEBOOK_DIR=$NBDIR/invisible_cities
	echo IC_NOTEBOOK_DIR set to $IC_NOTEBOOK_DIR
    
    	echo setting PYTHONPATH
	export PYTHONPATH=$ICTDIR
	echo PYTHONPATH is set to $PYTHONPATH
	fi
}



## Main command dispatcher

case $COMMAND in
    python_version_env)      python_version_env ;;
    notebook_env)            notebook_env ;;
    setup_git)               setup_git ;;
    
    *) echo Unrecognized command: ${COMMAND}
       echo
       echo Usage:
       echo
       echo source $0 python_version_env X.Y
       echo source $0 notebook_env
       echo source $0 setup_git
       ;;
esac
