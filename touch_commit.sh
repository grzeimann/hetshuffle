#!/bin/bash

# names
setting_file=hetdex_shuffle/__init__.py
property=touchsettings
tmp_commit=svn-commit.tmp
EDITOR=nano

function show_help(){
cat << EOF
Usage:
    ${0##*/} [-hnt]

Change the property '$property' of file '$setting_file' to a random string and
commit.
If the file '$tmp_commit' is available in the current or parent directory and
'-t' option is not given, it will be used to provide the commit message and, if
successfull, it will be removed.

    -h
        print this help
    -n
        don't do the commit
    -t
        don't use the '$tmp_commit' file for the svn commit message
EOF
}

# command line parsing
docommit=true
svntmp=true

while getopts "hnt" opt
do
  case "$opt" in
    h) show_help
       exit 0 ;;
    n) docommit=false ;;
    t) svntmp=false ;;
    '?') show_help
         exit 1 ;;
  esac
done

# set the property

# generate random 32 character alphanumeric string (upper and lowercase)
# from https://gist.github.com/earthgecko/3089509
NEW_UUID=$(LC_CTYPE=C tr -dc 'a-zA-Z0-9' < /dev/urandom| fold -w 32 | head -n 1)
svn propset $property $NEW_UUID $setting_file

# don't do the commit and exit
if [ "$docommit" = "false" ]
then
  echo "Skipping the commit as required"
  exit
fi

# look for the temporary svn commit file
if [ "$svntmp" = "false" ]
then
  tmp_commit=''
elif [ -f $tmp_commit ]
then
  tmp_commit=$tmp_commit
elif [ -f ../$tmp_commit ]
then
  tmp_commit=../$tmp_commit
else
  tmp_commit=''
fi

# do the commit
if [ -z "$tmp_commit" ]  # no temporary file
then
  svn commit
else
  echo Doing the commit with the old message
  cat "$tmp_commit"
  echo
  svn commit -F "$tmp_commit"
  if [ $? -eq 0 ]
  then
    rm "$tmp_commit"
    echo Temporary commit message "$tmp_commit" removed
  fi
fi
