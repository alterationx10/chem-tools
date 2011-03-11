#!/bin/bash

# Set your fortran compiler here
FC=gfortran

# This should be run in the main directory
WORK_DIR=`pwd`
BIN_DIR=$WORK_DIR/bin

function build-tools {
	cd $WORK_DIR/$1
	# Take care of any *.f files
	for F_FILES in `ls *.f 2> /dev/null`
	do
		$FC $F_FILES -o $BIN_DIR/`basename $F_FILES .f`
	done
	# Take care of any *.f90 files
	for F_FILES in `ls *.f90 2> /dev/null`
	do
		$FC $F_FILES -o $BIN_DIR/`basename $F_FILES .f90`
	done
	cd $WORK_DIR
}

while [ "$1" != "" ];
do
	# If you like tab complete, then you'll love this ;-)
	TARGET=`echo $1 | sed 's/[/]//g'`
	case $TARGET in 
		-h)		echo " "
				echo "Call with ./build-tools.sh target [target 2 ...]"
				echo "Target is a directory of tools you want to build"
				echo "i.e. ./build-tools.sh xyz-tools"
				echo " "
				echo "use ./build-tools.sh all to build all tools"
				echo " "
				echo "use ./build-tools.sh clean to delete all compiled tools"
				echo " "
				exit 0;;

		xyz-tools)	build-tools xyz-tools
				shift 1;;

		spectra-tools)	build-tools spectra-tools
				shift 1;;

		line-shape-generators)	
				build-tools line-shape-generators
				shift 1;;

		misc)		build-tools misc
				shift 1;;

		all)		build-tools xyz-tools
				build-tools spectra-tools
				build-tools line-shape-generators
				build-tools misc
				exit 0;;

		clean)		if [ "$BIN_DIR" != "" ]; then
					rm -f $BIN_DIR/*
				fi
				exit 0;;

		*)		echo "Need help? Run with -h"
				exit 0;;
	esac
done
