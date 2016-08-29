#!/bin/bash

if [ -e 'default.config' ]
then
    echo "Configs have been copied already"
else
    cp $baseDir/*.config ./
fi


