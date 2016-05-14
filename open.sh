#!/bin/bash

if type firefox &> /dev/null;
then
  firefox $1;
else
  open $1;
fi
