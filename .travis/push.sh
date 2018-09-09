#!/bin//sh

git_config() {
    git config --global user.email "mdsherman@betteridiot.tech"
    git config --global user.name "Marcus D. Sherman"
}

git_push() {
    git remote -v
    echo "Adding Remote"
    git remote add newOrigin https://${GH_TOKEN}@github.com/betteridiot/bamnostic.git
    echo "Fetching origin"
    git fetch newOrigin
    echo "Branches"
    git branch -vv
    echo "Checking out master"
    git checkout -b om newOrigin/master
    echo "Attempting merge"
    git merge newOrigin/devel
    echo "Attempting push"
    git push om
}

git_config
git_push