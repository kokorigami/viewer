openBrowser() {
  if which firefox &>/dev/null
  then
    firefox $1
  else
    open $1
  fi
}

openBrowser $1
