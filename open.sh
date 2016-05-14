function openBrowser {
  if which firefox
  then
    firefox $1
  else
    open $1
  fi
}

openBrowser $1
