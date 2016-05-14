#!/bin/bash

if type firefox &> /dev/null;
then
  export DISPLAY=:99.0;
  firefox $1;
else
  open $1;
fi
