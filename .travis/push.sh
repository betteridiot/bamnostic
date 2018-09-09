#!/bin//sh

git_config() {
    git config --global user.email "mdsherman@betteridiot.tech"
    git config --global user.name "Marcus D. Sherman"
}

git_push() {
    git remote add newOrigin https://${GH_TOKEN}@github.com/betteridiot/bamnostic.git
    git fetch newOrigin
    git pull newOrigin master
    git push newOrigin HEAD:master
}

git_config
git_push