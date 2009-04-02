#!/bin/bash

DIANA_HOME=/home/gdorea

add_to_path(){
	  if [ -d $1 ] ; then
      if [ $PATH ] ; then
	PATH=$1:"${PATH}"
      else
	PATH=$1
      fi
    fi
	export PATH
}

add_to_library_path () {
    if [ -d $1 ] ; then
      if [ $LD_LIBRARY_PATH ] ; then
	LD_LIBRARY_PATH=$1:"${LD_LIBRARY_PATH}"
      else
	LD_LIBRARY_PATH=$1
      fi
    fi
    export LD_LIBRARY_PATH
}

ompi_info=${DIANA_HOME}/ompi/bin/ompi_info

if [ -f $ompi_info ]; then
	OMPIBIN=`$ompi_info -path     bindir  -parsable | cut -d: -f3`
	OMPILIB=`$ompi_info -path     libdir  -parsable | cut -d: -f3`
	add_to_path ${OMPIBIN}
	add_to_library_path ${OMPILIB}
	unset  ompi_info OMPIBIN OMPILIB
fi

octave_config=${DIANA_HOME}/octave/bin/octave-config
if [ -f $octave_config ]; then
	OCTBIN=`$octave_config --print BINDIR`
	OCTAVE=`$octave_config --print PREFIX`
	add_to_path ${OCTBIN}
	export MPITB_HOME=${OCTAVE}/mpitb
fi

alias octave='octave -q'

mpirun octave --eval master