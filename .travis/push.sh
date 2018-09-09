#!/bin//sh

git_config() {
    git config --global user.email "mdsherman@betteridiot.tech"
    git config --global user.name "Marcus D. Sherman"
}

git_push() {
    echo "Adding Remote"
    git remote add origin https://${GH_TOKEN}@github.com/betteridiot/bamnostic.git > /dev/null 2>&1
    echo "Fetching origin"
    git fetch origin
    echo "Checking out master"
    git checkout -b master origin/master
    echo "Attempting merge"
    git merge origin/devel
    echo "Attempting push"
    git push origin/master
}

git_config
git_push