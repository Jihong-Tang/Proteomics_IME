#!/bin/bash

TIMESTAMP=$(date)
sed -i '' "/^Last updated:/s/.*/Last updated: $TIMESTAMP/" README.md