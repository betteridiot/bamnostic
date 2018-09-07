#!/bin//sh

git_config() {
    git config --global user.email "mdsherman@betteridiot.tech"
    git config --global user.name "Marcus D. Sherman"
}

git_push() {
    echo "Adding Remote"
    git remote add origin https://${GH_TOKEN}@github.com/betteridiot/bamnostic.git > /dev/null 2>&1
    echo "Fetching remote branch: devel"
    git branch --track dev origin/devel
    git branch
    git merge devel
    echo "Successful Merge"
    git push --quiet --set-upstream origin master
    echo "Push Complete"
}

git_config
git_push