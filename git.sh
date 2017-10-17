#!/bin/bash

read -p "Add, Commit, Push remote? "

if [[ $REPLY =~ ^[Yy]$ ]]
then
	echo    # (optional) move to a new line

    read -p "Provide commit message: "
	# Set git to use the credential memory cache
	#git config --global credential.helper cache

	# Set the cache to timeout after 4 hours (setting is in seconds)
	#git config --global credential.helper 'cache --timeout=14400'

	#stage new and modified files
	git add .

	#Record changes to the repository with the commit message
	git commit -m "$REPLY"

	#Push local changes to the github repository
	git push origin master
fi



