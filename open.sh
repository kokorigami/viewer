openBrowser() {
  if which firefox &>/dev/null
  then
    export DISPLAY=:0
    firefox $1
  else
    open $1
  fi
}

openBrowser $1
