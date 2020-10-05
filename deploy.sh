#!/bin/bash

scriptpath="$( cd "$(dirname "$0")" ; pwd -P )"

# make sure we are in the repo root
cwd=`pwd`

if [ $cwd != $scriptpath ]
then
    echo "Error: deploy.sh must be run from the repository root."
    exit 1
fi

# make sure we are on main branch
branch=`git rev-parse --abbrev-ref HEAD`

if [ "$branch" != "main" ]
then
    echo "Error: Cannot deploy from branch '$branch'. Switch to 'main' before deploying."
    exit 1
fi

# make sure there are no changes to commit
if git diff-index --quiet HEAD --
then
    msg=`git log -1 --pretty=%B`

    git pull origin main

    # build the site
    cd _site
    git fetch --quiet origin
    git reset --quiet --hard origin/main
    cd ..
    if ! bundle exec jekyll build --trace; then
        echo "Jekyll build failed. Master not updated."
        exit 1
    fi
    cd _site

    untracked=`git ls-files --other --exclude-standard --directory`

    # check if there are any changes on master
    if git diff --exit-code > /dev/null && [ "$untracked" = "" ]
    then
        echo "Nothing to update on master."
        cd ..
    else
        # deploy the static site
        git add . && \
        git commit -am "$msg" && \
        git push --quiet origin HEAD:master
        echo "Successfully built and pushed to master."
        cd ..
    fi

    # deploy main
    git push --quiet origin main
    echo "Deployment complete."
else
    echo "Error: Uncommitted main changes. Please commit or stash before updating master."
    exit 1
fi
